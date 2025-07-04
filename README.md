# DOC-WAN
Detecting Overlapping Communities in Weighted Attributed Networks
To run this algorithm, we need two files one containing the network in edgelist format, and the second one the attributes of the nodes.
First the code needs to be compiled using the following command in the directory where the code of DOC-WAN lies. 

g++ -std=c++11 doc-wan.cpp -o docwan

Once the executable file docwan is created, it can be used as follows:

./docwan ./network.txt ./node-attributes.txt -w -ov 0.5 -rh 0.5

Here network.txt is the file containing the network in edgelist format, and node-attributes.txt is the file containing the node attributes of the network. In the code all the attributes are assumed to be of floating point type. Apart from the two files, the code takes three more parameters, as shown above by three flags, -w, -ov, -rh. The flag -w indicates that the network is weighted. If the network is unweighted, it can be ommited. The flag -ov indicates the threshold value of the overlap, meaning that any two communities having overlap higher than or equal to this value are merged into one community. Likewise, -rh indicates the threshold value of the neighbourhood proximity. If the user is not sure about the values of these parameters, he or she can simply call the code for weighted networks as follows:

./docwan ./network.txt ./node-attributes.txt -w            

If the network is unweighted, then the code can be called as follows:

./docwan ./network.txt ./node-attributes.txt 

The output is written in the file called doc-wan-coms.txt, which resides in the same directory where the code lies.
