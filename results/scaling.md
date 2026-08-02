### Point lookup, indexed versus not

| n | Unindexed | Indexed | Speed-up |
| ---: | ---: | ---: | ---: |
| 1,000 | 37.4 µs | 7.2 µs | 5.2× |
| 4,000 | 128.5 µs | 7.1 µs | 18.1× |
| 16,000 | 496.1 µs | 7.5 µs | 66.2× |
| 64,000 | 1,986.1 µs | 7.3 µs | 271.9× |
| 256,000 | 9,649.0 µs | 8.2 µs | 1180.7× |

### Interval overlap, three structures

| n | B-tree scan | UCSC binning | R*Tree | Fastest |
| ---: | ---: | ---: | ---: | --- |
| 1,000 | 26.9 µs | 18.7 µs | 10.3 µs | rtree |
| 4,000 | 81.0 µs | 30.0 µs | 15.5 µs | rtree |
| 16,000 | 352.7 µs | 74.7 µs | 32.8 µs | rtree |
| 64,000 | 1,359.8 µs | 288.4 µs | 106.3 µs | rtree |
| 256,000 | 6,171.0 µs | 2,406.8 µs | 443.3 µs | rtree |
