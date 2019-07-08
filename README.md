# HOGWARTs


<a name="intro"/>

## Introduction

HOGWARTs (Hunt of Gravitational Wave Areas for Rapid Transients) is an algorithm which can be used to generate ranked lists of candidate galaxies within LIGO/Virgo localisation regions. The algorithm performs the following operations :

* Listen for a GCN containing a new gravitational wave alert
* Download and read the skymap using healpy
* Identify the contours containing 50%, 90% and 99% of the probability in the map
* Identify the galaxies within the GLADE V2 galaxy catalogue which lie within these localisation regions by querying Vizier
* Calculate a probability of the galaxy being associated with the gravitational wave source using the prioritisation algorithm outlined by [Arcavi et al. 2017](https://arxiv.org/abs/1710.05842):
	* The location probability measure is given as :

		<img src="https://latex.codecogs.com/gif.latex?S_{loc}=p_{loc}\,p_{dist}" title="S_{loc}=p_{loc}\,p_{dist}" />

	The probability that the GW source is at a certain location, <img src="https://latex.codecogs.com/gif.latex?p_{loc}" title="p_{loc}" />, is obtained from the pixel at the position of the galaxy in the sky map. 

	The distance to the merger computed by the LVC is contained in the pixel at the position of the galaxy. It is compared to the distance of the galaxy extracted from the filtered GLADE V2 galaxy catalogue to calculate the distance probability measure <img src="https://latex.codecogs.com/gif.latex?p_{dist}" title="p_{dist}" />:

	<img src="https://latex.codecogs.com/gif.latex?p_{&space;dist}=N_{&space;dist}\,exp(\frac{-[D-\mu_{&space;dist}]^2}{2\sigma_{&space;dist}^2})" title="p_{ dist}=N_{ dist}\,exp(\frac{-[D-\mu_{ dist}]^2}{2\sigma_{ dist}^2})" />

	where <img src="https://latex.codecogs.com/gif.latex?N_{dist}" title="N_{dist}" /> is a normalising factor, <img src="https://latex.codecogs.com/gif.latex?/=\mu_{dist}" title="\mu_{dist}" />  is the distance estimate and <img src="https://latex.codecogs.com/gif.latex?\sigma_{dist}" title="\sigma_{dist}" /> is the distance error computed by BAYESTAR and contained in the pixel at the galaxy's position sky map. D is distance to the galaxy from the filtered GLADE V2 catalogue. 

	* Short GRBs are found in the most massive galaxies, and B luminosity is a proxy for galaxy mass. Brighter galaxies are assigned a larger probability. The B-band luminosity is calculated using the apparent B magnitude and distance from the filtered GLADE V2 catalogue. This is used to calculate the Luminosity probability measure <img src="https://latex.codecogs.com/gif.latex?S_{lum}" title="S_{lum}" /> :

		<img src="https://latex.codecogs.com/gif.latex?S_{lum}=\frac{L_{B}}{\sum&space;L_{B}}" title="S_{lum}=\frac{L_{B}}{\sum L_{B}}" />

	* The overall probability of the merger occurring in a galaxy is given by 
	<img src="https://latex.codecogs.com/gif.latex?S=S_{loc}\,S_{lum}" title="S=S_{loc}\,S_{lum}" />
	This score computed for all galaxies and then the scores are normalised to add to 1.

* The galaxies are ranked in descending order and their properties are saved in a json and ascii file for each of the 50%, 90% and 99% localisation regions.
