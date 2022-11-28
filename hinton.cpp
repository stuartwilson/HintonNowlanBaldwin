#include <morph/HdfData.h>
#include <morph/Config.h>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
using namespace std;

morph::RandUniform<double, std::mt19937>* rng;

vector<int> randPermute(int x){
    vector<int> X(x);
    for(int i=0;i<x;i++){
        X[i] = i;
    }
    vector<int> Y(x);
    for(int i=0;i<x;i++){
        int index = floor(rng->get()*X.size());
        Y[i] = X[index];
        X.erase (X.begin()+index);
    }
    return Y;
}

std::vector<int> argSort(std::vector<double> x){
        // returns indices into x (descending)
        std::vector<int> sortedID(0);
        while(x.size()){
                float maxVal = -1e9;
                int maxID = -1;
                for(int i=0;i<x.size();i++){
                        if(x[i]>=maxVal){
                                maxVal = x[i];
                                maxID = i;
                        }
                }
                sortedID.push_back(maxID);
                x.erase(x.begin()+maxID);
        }
        return sortedID;
}

// ORIGINAL HINTON AND NOWLAN FITNESS FUNCTION
double getFitA(std::vector<int> x, vector<int> target, double f_min, int T){

    int N = x.size();
    for(int i=0;i<N;i++){
        if((x[i]!=target[i]) && (x[i]!=2)){
            return 1.0;
        }
    }
    int c = 0;
    for(int i=0;i<N;i++){
        if(x[i]==2){ c++; }
    }
    if(c){
        for(int t=0;t<T;t++){
	    if (rng->get()<pow(2,-(double)c)){
                return 1.0+(double)(N-1)*(T*1.-t)/(T*1.);
            }
        }
        return 1.0;
    } else {
        return (double)N;
    }
}


// SIMPLIFIED 'BINARY FITNESS' FUNCTION
double getFitB(std::vector<int> x, vector<int> target, double f_min, int T){

    int N = x.size();
    for(int i=0;i<N;i++){
        if((x[i]!=target[i]) && (x[i]!=2)){
            return f_min;
        }
    }
    int c = 0;
    for(int i=0;i<N;i++){
        if(x[i]==2){ c++; }
    }
    	double P = 1.0-pow((1.0-pow(0.5,(double)c)),(double)T);
    	if(rng->get()<P){
		return 1.0;
    	} else {
		return f_min;
    	}
}

// CONTINUOUS 'HYBRID' FITNESS FUNCTION
double getFitC(std::vector<int> x, vector<int> target, int T){

        int N = x.size();
        int undecided = 0;
        int correct = 0;
        for(int i=0;i<N;i++){
           if(x[i]==2){
                undecided++;
           } else if(x[i]==target[i]){
                correct++;
           }
        }
        double Q = 1.0-pow((1.0-pow(0.5,(double)undecided)),(double)T);
        return (correct*1.0+undecided*Q)/(double)N;

}

// FIXED 'LEARNING COST' FITNESS FUNCTION
double getFitD(std::vector<int> x, vector<int> target, double penaltyFactor){

	int N = x.size();
	int undecided = 0;
	int correct = 0;
	for(int i=0;i<N;i++){
	   if(x[i]==2){
		undecided++;
	   } else if(x[i]==target[i]){
		correct++;
	   }
	}
	return (correct*1.0+undecided*penaltyFactor)/(double)N;
}

// CHUNK VERSION OF BINARY FITNESS
double getFitE(std::vector<int> x, vector<int> target, double f_min, int _T, int chunks){

    int T = _T/chunks;
    int N = x.size();
    int nPerChunk = N/chunks;
    vector<double> Fs(chunks,0.0);
    for(int k=0;k<chunks;k++){
	int startOfChunk = chunks*nPerChunk;
	for(int l=0;l<nPerChunk;l++){
		int i=startOfChunk+l;
		 if((x[i]!=target[i]) && (x[i]!=2)){
            		Fs[k]=f_min;
			break;
        	 }
	}
       	int c = 0;
    	for(int l=0;l<nPerChunk;l++){
		int i=startOfChunk+l;
        	if(x[i]==2){ 
			c++; 
		}
    	}

        double P = 1.0-pow((1.0-pow(0.5,(double)c)),(double)T);
        if(rng->get()<P){
                Fs[k] = 1.0;
		break;
        } else {
                Fs[k] = f_min;
		break;
        }
    }
    float f = 0.0;
    for(int k=0;k<chunks;k++){
	f+=Fs[k];
    }
    return f/(float)chunks;
}

// HINTON AND NOWLAN SELECTION METHOD (SELECTION PROBABILITY BASED ON RELATIVE FITNESS)
vector<vector<int> > selectionHN (vector<vector<int> > X, vector<double> F, int crossOverPoint){

	int P = X.size();
	int N = X[0].size();

        vector<double> lower (P,0.);
        vector<double> upper (P,0.);

        double cumSum = 0.0;
        for(int p=0;p<P;p++){
                lower[p] = cumSum;
                cumSum += F[p];
                upper[p] = cumSum;
        }

        vector<vector<int> > Mum(P,vector<int>(N,0));
        for(int p=0;p<P;p++){
            double r = rng->get()*cumSum;
            for(int q=0;q<P;q++){
                if ((lower[q]<=r) && (r<upper[q])){
                    Mum[p] = X[q];
                    break;
                }
            }
        }

        vector<vector<int> > Dad(P,vector<int>(N,0));
        for(int p=0;p<P;p++){
            double r = rng->get()*cumSum;
            for(int q=0;q<P;q++){
                if ((lower[q]<=r) && (r<upper[q])){
                    Dad[p] = X[q];
                    break;
                }
            }
        }

        X = Mum;


	if(crossOverPoint==0){
	/*
   		for(int p=0;p<P;p++){
            		int crossOver = rng->get()*N;
            		for(int i=crossOver;i<N;i++){
                		X[p][i] = Dad[p][i];
            		}
        	}
        */
         for(int p=0;p<P;p++){
            		int crossOver = rng->get()*(N-1);
            		for(int i=crossOver+1;i<N;i++){
                		X[p][i] = Dad[p][i];
            		}
        	}

	} else {
        	for(int p=0;p<P;p++){
            		for(int i=0;i<crossOverPoint;i++){
                		X[p][i] = Dad[p][i];
            		}
        	}
	}

	return X;
}

// MACKAY SELECTION METHOD (PAIRS FROM TOP HALF EACH SEED FOUR OFFSPRING)
vector<vector<int> >  selectionMK (vector<vector<int> > X, vector<double> F, int crossOverPoint){

	int P = X.size();
	int N = X[0].size();
	int nHalf = (int)(P/2);
        vector<vector<int> > Xcp = X;
	std::vector<int> shuffleIndex = randPermute(P);
	std::vector<double> shuffleVal(P);
	for(int i=0;i<P;i++){
		shuffleVal[i]=F[shuffleIndex[i]];
	}
        std::vector<int> sortedID = argSort(shuffleVal);
        std::vector<int> permutation = randPermute(nHalf);
        int k=0;
        for(int i=0;i<nHalf;i+=2){
                int mum = shuffleIndex[sortedID[permutation[i]]];
                int dad = shuffleIndex[sortedID[permutation[i+1]]];
                for(int j=0;j<4;j++){
                        X[k] = Xcp[mum];
			if(crossOverPoint==0){
                        	for(int l=0;l<N;l++){
                                	if(rng->get()<0.5){
                                        	X[k][l] = Xcp[dad][l];
                                	}
                        	}
			} else {
 				for(int l=0;l<crossOverPoint;l++){
                                      	X[k][l] = Xcp[dad][l];
                                }
			}
                        k++;
                }
        }
	return X;
}



int main (int argc, char **argv){

//	std::cout<<"A"<<std::endl;
    if (argc < 3) { std::cerr << "\nUsage: ./hintonBatchVolatileConfigurable configfile logdir genomeMutationRate targetMutationRate numMutable seed(optional)"; return -1; }

// std::cout<<"B"<<std::endl;
    srand(time(NULL)); // note may not be different for simultaneously launched progs
    int seed = rand();
    if(argc==7){
        seed = std::stoi(argv[6]);
    }
    morph::RandUniform<double, std::mt19937> _rng(seed);
    rng = &_rng;

//std::cout<<"C"<<std::endl;
    std::string paramsfile (argv[1]);
    morph::Config conf(paramsfile);
    if (!conf.ready) { std::cerr << "Error setting up JSON config: " << conf.emsg << std::endl; return 1; }

//std::cout<<"D"<<std::endl;

    std::string logpath = argv[2];
    std::ofstream logfile;
    morph::Tools::createDir (logpath);
    { std::stringstream ss; ss << logpath << "/log.txt"; logfile.open(ss.str());}
    logfile<<"Hello."<<std::endl;


//std::cout<<"E"<<std::endl;
    int selectionMode = conf.getInt("selectionMode", 0);
    int fitnessMode = conf.getInt("fitnessMode", 0);

    int G = conf.getInt("generations", 1000);
    int P = conf.getInt("populationSize", 1000);
    int T = conf.getInt("lifetimeLearningTrials", 1000); // NOTE SET T=-1 IN CONFIG TO USE PARAM LINE INPUT ARGV[4] TO SPECIFY IT INSTEAD    
    int N = conf.getInt("numAlleles", 20);
    double initialQs = conf.getDouble("initialProportionLearnable", 0.5);
    int crossOverPoint = conf.getInt("crossOverPoint", 0);

    int mutationRateRes = conf.getInt("genomeMutationRateRes", 50);
    int targetMutationRateRes = conf.getInt("targetMutationRateRes", 50);

    int chunks = conf.getInt("chunks",1);

    double mutationRate = (double)std::stoi(argv[3])/(double)(mutationRateRes-1);
    
    double targetMutationRate; 
    if(T<0){
    	T = stoi(argv[4]);
    	targetMutationRate = 0.0;
    } else {
    	targetMutationRate = (double)std::stoi(argv[4])/(double)(targetMutationRateRes-1);
    }
    
    int numMutable = std::stoi(argv[5]);
    double learningPenaltyFactor = conf.getDouble("learningPenaltyFactor", 0.5);
    double norm = 1.0/((double)P*(double)N);

    vector<int> rp= randPermute(N);
    vector<int> mutableID;
    for(int i=0;i<numMutable;i++){
	mutableID.push_back(rp[i]);
    }

    double f_min = conf.getDouble("f_min", 1.0/(double)P);

    vector<int> target(N,1);

    // setup initial population
    vector<vector<int> > X(P,vector<int>(N,0));
    for(int p=0;p<P;p++){
        for(int i=0;i<N;i++){
            if (rng->get()<initialQs){
                X[p][i]=2;
            }
            else if (rng->get()<0.5){
                X[p][i]=1;
            }
        }
    }

    vector<double> p0, p1, p2, pCorr;
 
    // MAIN EVOLUTINOARY LOOP

    for(int g=0;g<G;g++){

        // extract states
        vector<int> counts(3,0);
        int corr = 0;
        for(int p=0;p<P;p++){
            for(int i=0;i<N;i++){
                counts[X[p][i]]++;
            }
            for(int i=0;i<N;i++){
                if(target[i]==X[p][i]){
                    corr++;
                }
            }
        }
        p0.push_back(counts[0]*norm);
        p1.push_back(counts[1]*norm);
        p2.push_back(counts[2]*norm);
        pCorr.push_back((double)corr*norm);


	// vary target
	 for(int i=0;i<numMutable;i++){
            if(rng->get()<targetMutationRate){
                target[mutableID[i]] = 1-target[mutableID[i]];
	    }
        }

	// compute fitness
	vector<double> F(P,1.);
	switch(fitnessMode){
		case(0): {
                        for(int p=0;p<P;p++){
                                F[p] = getFitA(X[p],target,f_min,T);
                        }
                } break;
   		case(1): {
                         for(int p=0;p<P;p++){
                                F[p] = getFitB(X[p],target,f_min,T);
                        }
                } break;
		case(2): {
        		for(int p=0;p<P;p++){
            			F[p] = getFitC(X[p],target,T);
        		}
		} break;
		case(3): {
			 for(int p=0;p<P;p++){
                                F[p] = getFitD(X[p],target,learningPenaltyFactor);
                        }
		} break;
		 case(4): {
			vector<int> x(N/chunks,0);
			vector<int> t(N/chunks,0);
			for(int p=0;p<P;p++){
				double f=0.0;
				int k=0;
				for(int i=0;i<chunks;i++){
					for(int j=0;j<N/chunks;j++){
						x[j] = X[p][k];
						t[j] = target[k];
						k++;
					}
					f+= getFitB(x,t,f_min,T/chunks);
				}
				F[p] = f/(float)chunks;
			}
                         //for(int p=0;p<P;p++){
                         //       F[p] = getFitE(X[p],target,f_min,T,chunks);
                        //}
                } break;
	}


	// selection
	switch(selectionMode){
		case(0): {
    			X = selectionHN(X,F,crossOverPoint);
		} break;
		case(1): {
			X = selectionMK(X,F,crossOverPoint);
		} break;
	}

	// mutation
	if(mutationRate>0){
            for(int p=0;p<P;p++){
                for(int i=0;i<N;i++){
                    if(rng->get()<mutationRate){
                        X[p][i] = floor(rng->get()*3);
                    }
                }
            }
        }



    }


    // save states to file
    std::stringstream fname;
    fname << logpath << "/out.h5";
    morph::HdfData data(fname.str());
    std::stringstream path;

    path.str(""); path.clear(); path << "/p0";
    data.add_contained_vals (path.str().c_str(), p0);

    path.str(""); path.clear(); path << "/p1";
    data.add_contained_vals (path.str().c_str(), p1);

    path.str(""); path.clear(); path << "/p2";
    data.add_contained_vals (path.str().c_str(), p2);

    path.str(""); path.clear(); path << "/pCorr";
    data.add_contained_vals (path.str().c_str(), pCorr);

    return 0;
}


