### Point lookup, indexed versus not

| n | Unindexed | Indexed | Speed-up |
| ---: | ---: | ---: | ---: |
| 1,000 | 30.4 µs | 7.9 µs | 3.8× |
| 4,000 | 105.6 µs | 7.0 µs | 15.0× |
| 16,000 | 389.4 µs | 7.5 µs | 51.6× |
| 64,000 | 1,520.4 µs | 7.5 µs | 203.7× |
| 78,733 | 2,336.0 µs | 8.7 µs | 270.0× |

### Interval overlap, three structures

| n | B-tree scan | UCSC binning | R*Tree | Fastest |
| ---: | ---: | ---: | ---: | --- |
| 1,000 | 30.3 µs | 16.4 µs | 8.8 µs | rtree |
| 4,000 | 95.8 µs | 17.6 µs | 9.8 µs | rtree |
| 16,000 | 364.7 µs | 24.1 µs | 11.6 µs | rtree |
| 64,000 | 1,383.4 µs | 36.6 µs | 20.1 µs | rtree |
| 78,733 | 1,615.1 µs | 40.6 µs | 22.0 µs | rtree |
