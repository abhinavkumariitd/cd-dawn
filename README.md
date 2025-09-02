#CD-DAWN
Community Detection in Directed And Weighted Networks (CD-DAWN)
To run this algorithm, a network file containing the directed network in edge list format is required.
First, the code must be compiled on a Linux/Unix machine using the following command in the directory where the code of CD-DAWN lies. 

g++ -std=c++11 cd-dawn.cpp -o cd-dawn

Once the executable file cd-dawn is created, it can be used as follows:

./cd-dawn ./network.txt -w -rh [option] -ov [option]

The algorithm takes a network file, network.txt, as input, representing a directed network in edge list format—either with two columns for unweighted or three for weighted networks. It also accepts three flags: -w, -ov, and -rh. The flag -w indicates that the network is weighted; if the network is unweighted, this flag can be omitted. The flag -ov specifies the overlap threshold, meaning that any two communities with an overlap greater than or equal to this value are merged into one community (default: 0.40). Likewise, -rh indicates the threshold value for neighbourhood proximity (default:0.30).


If the user is not sure about the values of these parameters, he or she can simply call the code for weighted networks as follows:

./cd-dawn ./network.txt  -w            

If the network is unweighted, then the code can be called as follows:

./cd-dawn ./network.txt 

The output is written in the file called cd-dawn-coms.txt, which resides in the same directory where the code lies.

If you use this code for research purpose, please cite the following article:

@article{kumar2025overlapping,
  title={Overlapping community detection with a new modularity measure in directed weighted networks},
  author={Kumar, Abhinav and Kumari, Anjali and Kumar, Pawan and Dohare, Ravins},
  journal={Data Mining and Knowledge Discovery},
  volume={39},
  number={6},
  pages={1--37},
  year={2025},
  publisher={Springer}
}
