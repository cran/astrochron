\name{getLaskar}
\alias{getLaskar}
\title{Download Laskar et al. (2004, 2011a, 2011b) astronomical solutions}
\description{
Download Laskar et al. (2004, 2011a, 2011b) astronomical solutions.
}
\usage{
getLaskar(sol="la04")
}

\arguments{
  \item{sol}{A character string that specifies the astronomical solution to download: "la04","la10a","la10b","la10c","la10d","la11"}
}

\details{
la04 : three columns containing precession angle, obliquity, and eccentricity of Laskar et al. (2004)

la10a : one column containing the la10a eccentricity solution of Laskar et al. (2011a)

la10b : one column containing the la10b eccentricity solution of Laskar et al. (2011a)

la10c : one column containing the la10c eccentricity solution of Laskar et al. (2011a)

la10d : one column containing the la10d eccentricity solution of Laskar et al. (2011a)

la11 : one column containing the la11 eccentricity solution of Laskar et al. (2011b; please also cite 2011a)
}

\references{
J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B. Levrard, B., 2004, 
\emph{A long term numerical solution for the insolation quantities of the Earth}: 
Astron. Astrophys., Volume 428, 261-285.

Laskar, J., Fienga, A., Gastineau, M., Manche, H., 2011a, 
\emph{La2010: A new orbital solution for the long-term motion of the Earth}: Astron. Astrophys., Volume 532, A89.

Laskar, J., Gastineau, M., Delisle, J.-B., Farres, A., Fienga, A.: 2011b,
\emph{Strong chaos induced by close encounters with Ceres and Vesta, Astron}: Astrophys., Volume 532, L4.
}