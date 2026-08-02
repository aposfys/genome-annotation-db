### Point lookup, indexed versus not

| n | Unindexed | Indexed | Speed-up |
| ---: | ---: | ---: | ---: |
| 1,000 | 19.8 µs | 3.8 µs | 5.2× |
| 4,000 | 68.1 µs | 4.3 µs | 15.9× |
| 16,000 | 260.7 µs | 3.7 µs | 69.6× |
| 64,000 | 1,040.4 µs | 3.8 µs | 275.8× |
| 256,000 | 4,930.8 µs | 3.8 µs | 1284.2× |

### Interval overlap, three structures

| n | B-tree scan | UCSC binning | R*Tree | Fastest |
| ---: | ---: | ---: | ---: | --- |
| 1,000 | 15.0 µs | 10.0 µs | 5.4 µs | rtree |
| 4,000 | 44.7 µs | 15.9 µs | 7.8 µs | rtree |
| 16,000 | 175.9 µs | 38.5 µs | 15.0 µs | rtree |
| 64,000 | 690.9 µs | 152.7 µs | 49.9 µs | rtree |
| 256,000 | 3,178.1 µs | 1,195.1 µs | 184.3 µs | rtree |
