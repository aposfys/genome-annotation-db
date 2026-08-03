### Point lookup, indexed versus not

| n | Unindexed | Indexed | Speed-up |
| ---: | ---: | ---: | ---: |
| 1,000 | 17.7 µs | 3.6 µs | 4.9× |
| 4,000 | 54.0 µs | 3.7 µs | 14.4× |
| 16,000 | 195.6 µs | 3.8 µs | 51.5× |
| 64,000 | 781.7 µs | 4.2 µs | 185.2× |
| 78,733 | 1,179.4 µs | 3.9 µs | 305.4× |

### Interval overlap, three structures

| n | B-tree scan | UCSC binning | R*Tree | Fastest |
| ---: | ---: | ---: | ---: | --- |
| 1,000 | 14.1 µs | 8.2 µs | 4.4 µs | rtree |
| 4,000 | 46.9 µs | 8.9 µs | 5.0 µs | rtree |
| 16,000 | 190.4 µs | 10.6 µs | 5.9 µs | rtree |
| 64,000 | 742.9 µs | 19.7 µs | 10.2 µs | rtree |
| 78,733 | 885.4 µs | 20.7 µs | 11.1 µs | rtree |
