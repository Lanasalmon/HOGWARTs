# HOGWARTs


<a name="intro"/>

## Introduction

HOGWARTs (Hunt of Gravitational Wave Areas for Rapid Transients) is an algorithm which can be used to generate ranked lists of candidate galaxies within LIGO/Virgo localisation regions. The algorithm performs the following operations :

* Listen for a GCN containing a new gravitational wave alert
* Download and read the skymap using healpy
* Identify the contours containing 50%,90% and 99% of the probability in the map
* Identify the galaxies within the GLADE V2 galaxy catalogue which lie within these localisation regions by querying Vizier
* Calculate a probability of the galaxy being associated with the gravitational wave source using the prioritisation algorithm outlined by Arcavi et al. 2017 :
* The location probability measure is given as
$$
\begin{equation}
S_{loc}=p_{loc}\,p_{dist}
\end{equation}
$$

The probability that the GW source is at a certain location, p_{loc}, is obtained from the pixel at the position of the galaxy in the sky map. 

The distance to the merger computed by the LVC is contained in the pixel at the position of the galaxy. It is compared to the distance of the galaxy extracted from the filtered GLADE V2 galaxy catalogue to calculate the distance probability measure $p_{dist}$:
\begin{equation} 
p_{ dist}=N_{ dist}\,exp(\frac{-[D-\mu_{ dist}]^2}{2\sigma_{ dist}^2})
\end{equation}
where $N_{dist}$ is a normalising factor, $\mu_{dist}$ is the distance estimate and $\sigma_{dist}$ is the distance error computed by BAYESTAR and contained in the pixel at the galaxy's position sky map. D is distance to the galaxy from the filtered GLADE V2 catalogue. 

* Short GRBs are found in the most massive galaxies, and B luminosity is a proxy for galaxy mass. Brighter galaxies are assigned a larger probability. The B-band luminosity is calculated using the apparent B magnitude and distance from the filtered GLADE V2 catalogue. This is used to calculate the Luminosity probability measure $S_{lum}$:

\begin{equation}
S_{lum}=\frac{L_{B}}{\sum L_{B}}
\end{equation}

* The overall probability of the merger occurring in a galaxy is given by 
\begin{equation}
S=S_{loc}\,S_{lum}
\end{equation}

This score computed for all galaxies and then the scores are normalised to add to 1.


