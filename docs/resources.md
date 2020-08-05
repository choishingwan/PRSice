# Introduction
PRSice relies on a number of open source projects to achieve the current performance. 
We also used algorithm found in other projects and translate them into C++ code for our own use. Below are number of projects we relies on

## Open source projects

| Project | Developer(s) | Description|
|:-|:-:|:-|
|[PLINK 2](https://www.cog-genomics.org/software) |Christopher Chang | Provide the backbone of the clumping algorithm and PRS calculation|
|[BGEN lib](https://bitbucket.org/gavinband/bgen) |Gavin Band| Provide API to handle BGEN files. Slight modification were made to accomodate PRSice's usage|
|[Eigen C++](https://eigen.tuxfamily.org)| Ga&euml;l Guennebaud and Beno&icirc;t Jacob and others | For all matrix algebra|
|[gzstream](http://www.cs.unc.edu/Research/compgeom/gzstream/)|  Deepak Bandyopadhyay and Lutz Kettner| For reading gz files|
| [fastglm](https://github.com/jaredhuling/fastglm) |  	Jared Huling, Douglas Bates, Dirk Eddelbuettel, Romain Francois and Yixuan Qiu | Basis of our glm class| 
| [RcppEigen](https://github.com/RcppCore/RcppEigen) |   	Douglas Bates, Dirk Eddelbuettel, Romain Francois, and Yixuan Qiu | Provide the fastlm algorithm |

