#include <morph/HdfData.h>
#include <morph/Config.h>

#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
using namespace std;

morph::RandUniform<double, std::mt19937>* rng;

double getFit(std::vector<int> x, int T, double f_min, vector<int> target){

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
    if(c){
        for(int t=0;t<T;t++){
            double f = 1.0;
            for(int i=0;i<c;i++){
                f *= (rng->get()<0.5);
            }
            if (f==1){
                //return 1.0-(T*1.-t)/(T*1.);
                return (T*1.-t)/(T*1.);
            }
        }
        return f_min;
    } else {
        return 1.0;
    }
}

int main (int argc, char **argv){

    if (argc < 3) { std::cerr << "\nUsage: ./evodevo configfile logdir seed(optional)"; return -1; }

    srand(time(NULL)); // note may not be different for simultaneously launched progs
    int seed = rand();
    if(argc==4){
        seed = std::stoi(argv[3]);
    }
    morph::RandUniform<double, std::mt19937> _rng(seed);
    rng = &_rng;

    std::string paramsfile (argv[1]);
    morph::Config conf(paramsfile);
    if (!conf.ready) { std::cerr << "Error setting up JSON config: " << conf.emsg << std::endl; return 1; }

    std::string logpath = argv[2];
    std::ofstream logfile;
    morph::Tools::createDir (logpath);
    { std::stringstream ss; ss << logpath << "/log.txt"; logfile.open(ss.str());}
    logfile<<"Hello."<<std::endl;

    int G = conf.getInt("generations", 1000);   // generations
    int P = conf.getInt("populationSize", 1000);
    int T = conf.getInt("lifetimeLearningTrials", 1000);   // nodes
    int N = conf.getInt("numAlleles", 20);   // learning steps
    double initialQs = conf.getDouble("initialProportionLearnable", 0.5);
    double mutationRate = conf.getDouble("genomeMutationRate", 0.0);
    double targetMutationRate = conf.getDouble("targetMutationRate", 0.0);
    int Nover2 = N/2;

    double f_min = 1./(double)P;

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

    vector<double> lower (P,0.);
    vector<double> upper (P,0.);

    vector<double> p0, p1, p2, pCorr;
    double norm = 1.0/((double)P*(double)N);

    // main loop
    vector<double> F(P,1.);

    for(int g=0;g<G;g++){

        for(int i=0;i<N;i++){
            if(rng->get()<targetMutationRate){
                target[i] = 1-target[i];
            }
        }

        double sumF = 0.;
        for(int p=0;p<P;p++){
            F[p] = getFit(X[p],T,f_min,target);
            sumF += F[p];
        }

        double cumSum = 0.;
        for(int p=0;p<P;p++){
            lower[p] = cumSum;
            cumSum += F[p];
            upper[p] = cumSum;
        }

        // store states
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

        // recombination
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

        for(int p=0;p<P;p++){
            for(int i=Nover2;i<N;i++){
                X[p][i] = Dad[p][i];
            }
        }

        if(mutationRate>0){
            for(int p=0;p<P;p++){
                for(int i=0;i<N;i++){
                    if(rng->get()<mutationRate){
                        X[p][i] = floor(rng->get()*3);
                    }
                }
            }
        }

        std::cout<<"  "<<g<<"     \r"<<std::flush;

    }

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


