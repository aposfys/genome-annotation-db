"""Database connection and dialect handling.

The same schema and the same queries run on SQLite and MySQL. SQLite is the
default because it needs no server, which is what makes this repository
runnable by anyone and testable in CI; MySQL is supported because that is what
the database was originally built on and what a shared deployment would use.

Portability is handled here rather than by carrying two copies of every
statement. Only two things actually differ in practice: the parameter
placeholder style, and a handful of DDL keywords.
"""

from __future__ import annotations

import os
import re
import sqlite3
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlparse

DEFAULT_SQLITE_PATH = Path("data/genome.sqlite")

# Connection details are read from the environment. They are never written into
# the source: an earlier version of this project hardcoded a live host, user and
# password, which is exactly the mistake that puts credentials into git history.
DB_URL_VAR = "GENOMEDB_URL"


@dataclass(frozen=True)
class Backend:
    """Which database engine a connection speaks."""

    name: str  # "sqlite" or "mysql"

    @property
    def placeholder(self) -> str:
        return "?" if self.name == "sqlite" else "%s"

    @property
    def is_sqlite(self) -> bool:
        return self.name == "sqlite"


SQLITE = Backend("sqlite")
MYSQL = Backend("mysql")

# Statements written with SQLite's placeholder are rewritten for MySQL. Doing it
# in one place means queries can be authored in a single style.
_PLACEHOLDER = re.compile(r"\?")


def adapt(sql: str, backend: Backend) -> str:
    """Rewrite a statement for the target backend."""
    if backend.is_sqlite:
        return sql
    return _PLACEHOLDER.sub("%s", sql)


class Connection:
    """A thin wrapper that hides the two drivers' differences."""

    def __init__(self, raw: Any, backend: Backend) -> None:
        self.raw = raw
        self.backend = backend

    def execute(self, sql: str, params: Sequence[Any] = ()) -> Any:
        cursor = self.raw.cursor()
        cursor.execute(adapt(sql, self.backend), tuple(params))
        return cursor

    def executemany(self, sql: str, rows: Sequence[Sequence[Any]]) -> int:
        cursor = self.raw.cursor()
        cursor.executemany(adapt(sql, self.backend), [tuple(r) for r in rows])
        count = cursor.rowcount
        cursor.close()
        return count

    def executescript(self, sql: str) -> None:
        """Run a multi-statement script.

        Neither driver does this identically: sqlite3 has ``executescript``,
        while MySQL needs the statements fed one at a time.
        """
        if self.backend.is_sqlite:
            self.raw.executescript(sql)
            return

        cursor = self.raw.cursor()
        for statement in split_statements(sql):
            cursor.execute(statement)
        cursor.close()

    def query(self, sql: str, params: Sequence[Any] = ()) -> tuple[list[str], list[tuple]]:
        """Run a SELECT and return its column names and rows."""
        cursor = self.execute(sql, params)
        rows = cursor.fetchall()
        columns = [d[0] for d in cursor.description] if cursor.description else []
        cursor.close()
        return columns, [tuple(r) for r in rows]

    def scalar(self, sql: str, params: Sequence[Any] = ()) -> Any:
        _, rows = self.query(sql, params)
        return rows[0][0] if rows and rows[0] else None

    def count_vm_steps(self, sql: str, params: Sequence[Any] = ()) -> int:
        """Virtual-machine instructions the engine executes for a query.

        Wall-clock time answers "how fast is this machine"; this answers "how
        much work is this query", and gives the same number on every machine.
        A speed-up expressed in VM steps is therefore reproducible in a way a
        millisecond figure is not.

        Only SQLite exposes this cheaply, through the progress handler; on
        other backends it returns -1 rather than pretending.
        """
        if not self.backend.is_sqlite:
            return -1

        steps = 0

        def tick() -> int:
            nonlocal steps
            steps += 1
            return 0

        # A period of 1 counts every instruction. It is slow, which is why this
        # runs as its own pass and never inside a timed one.
        self.raw.set_progress_handler(tick, 1)
        try:
            cursor = self.raw.cursor()
            cursor.execute(adapt(sql, self.backend), tuple(params))
            cursor.fetchall()
            cursor.close()
        finally:
            self.raw.set_progress_handler(None, 0)
        return steps

    def commit(self) -> None:
        self.raw.commit()

    def close(self) -> None:
        self.raw.close()


def split_statements(script: str) -> list[str]:
    """Split a SQL script on semicolons, ignoring those inside strings or comments."""
    statements: list[str] = []
    buffer: list[str] = []
    in_single = in_double = in_line_comment = False

    for index, character in enumerate(script):
        if in_line_comment:
            buffer.append(character)
            if character == "\n":
                in_line_comment = False
            continue

        if not in_single and not in_double and script.startswith("--", index):
            in_line_comment = True
            buffer.append(character)
            continue

        if character == "'" and not in_double:
            in_single = not in_single
        elif character == '"' and not in_single:
            in_double = not in_double

        if character == ";" and not in_single and not in_double:
            statement = "".join(buffer).strip()
            if statement:
                statements.append(statement)
            buffer = []
            continue

        buffer.append(character)

    tail = "".join(buffer).strip()
    if tail:
        statements.append(tail)
    return statements


def connect(url: str | None = None, sqlite_path: Path | None = None) -> Connection:
    """Open a connection from a URL, the environment, or a local SQLite file.

    Args:
        url: ``mysql://user:password@host/database`` or ``sqlite:///path``. When
            omitted, ``$GENOMEDB_URL`` is used; failing that, a local SQLite file.
        sqlite_path: Where the default SQLite database lives.
    """
    url = url or os.environ.get(DB_URL_VAR)

    if url and url.startswith("mysql"):
        return _connect_mysql(url)

    path = Path(url[len("sqlite:///") :]) if url and url.startswith("sqlite:///") else None
    path = path or sqlite_path or DEFAULT_SQLITE_PATH
    return _connect_sqlite(path)


def _connect_sqlite(path: Path) -> Connection:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = sqlite3.connect(path)
    # Off by default in SQLite, and this schema depends on them.
    raw.execute("PRAGMA foreign_keys = ON")
    return Connection(raw, SQLITE)


def _connect_mysql(url: str) -> Connection:
    try:
        import mysql.connector
    except ImportError as error:  # pragma: no cover - optional dependency
        raise RuntimeError(
            "MySQL support needs mysql-connector-python: pip install '.[mysql]'"
        ) from error

    parts = urlparse(url)
    if not parts.hostname or not parts.path.lstrip("/"):
        raise ValueError(f"Malformed MySQL URL: {url}")

    raw = mysql.connector.connect(
        host=parts.hostname,
        port=parts.port or 3306,
        user=unquote(parts.username or ""),
        password=unquote(parts.password or ""),
        database=parts.path.lstrip("/"),
    )
    return Connection(raw, MYSQL)


@contextmanager
def connection(
    url: str | None = None, sqlite_path: Path | None = None
) -> Iterator[Connection]:
    """Open a connection and always close it."""
    conn = connect(url, sqlite_path)
    try:
        yield conn
    finally:
        conn.close()
