\name{tones}
\alias{tones}
\title{Calculate all possible difference and combinations tones}
\description{
Determine all possible difference and combinations tones from a set of frequencies, and
find the closest one to a specified frequency
}
\usage{
tones(a=NULL,freqs=NULL,f=T)
}
\arguments{
  \item{a}{The frequency you are seeking to match, in cycles/ka.}
  \item{freqs}{The vector of frequencies from which to calculate difference and combination tones, in cycles/ka.}
  \item{f}{Output results as frequencies (cycles/ka)? If false, will output results as periods (ka). (T or F)}
}