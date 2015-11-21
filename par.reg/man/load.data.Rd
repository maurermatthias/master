\name{load.data}
\alias{load.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Load Test-Dataset
%%  ~~function to do ... ~~
}
\description{ 
Loads a example dataset to test the functionality of the package, especially the function pr().
}
\usage{
load.data()
}
%- maybe also 'usage' for other objects documented here.
\details{ details test
%%  ~~ If necessary, more details than the description above ~~
}
\value{ 
A matrix with 200 rows and 3 columns is retrurned. Each row specifies an fatigue data observation, whereas the first column describes the stress-level (predictor), the second column describes the cyrcles to failure (observed value) and the third column describes if this obervation is a run-through observation.
}
\references{ 
The data is taken from C.-l. Shen, "The statistical analysis of fatigue data.," 1994, page 259.
}
\author{ Matthias Maurer
}
%\note{ note test
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{ seealso test
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{ 
##load data
data=load.data()

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
