\name{ma_deg}
\alias{ma_deg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MA pairwise plot generator using DESeq2
%%  ~~function to do ... ~~
}
\description{
This function takes in a list of pairwise comparison conditions from DESeq matrix data along with their annotations and outputs a list with the first entry being DESeq gene data and second entry containing a list of the MA-plot objects.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ma_deg(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{condits}{ A list of conditions from which to compare. Ex: "t0", "t1h_veh", "t1h_PDGF".
%%     ~~Describe \code{x} here~~
  }
  \item{annotations}{ The annotations of the conditions. Contains each condition in a list repeated n times, n equalling the number of repetitions in the study. Ex: n = 3, "t0","t0","t0","t1h_veh","t1h_veh","t1h_veh"
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
