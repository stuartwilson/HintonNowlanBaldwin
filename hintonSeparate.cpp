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


// ORIGINAL HINTON AND NOWLAN FITNESS FUNCTION
double getFitA(std::vector<std::vector<bool> > x, vector<bool> target, double f_min, int T){

    int N = x.size();
    for(int i=0;i<N;i++){
        if((x[i][0]!=target[i]) && (!x[i][1])){
            return 1.0;
        }
    }
    int c = 0;
    for(int i=0;i<N;i++){
        if(x[i][1]){ c++; }
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
double getFitB(std::vector<std::vector<bool> > x, vector<bool> target, double f_min, int T){

    int N = x.size();
    for(int i=0;i<N;i++){
        if((x[i][0]!=target[i]) && (!x[i][1])){
            return f_min;
        }
    }
    int c = 0;
    for(int i=0;i<N;i++){
        if(x[i][1]){ c++; }
    }
    double P = 1.0-pow((1.0-pow(0.5,(double)c)),(double)T);
    if(rng->get()<P){
        return 1.0;
    } else {
        return f_min;
    }
}


// HINTON AND NOWLAN SELECTION METHOD (SELECTION PROBABILITY BASED ON RELATIVE FITNESS)
vector<vector<vector<bool> > > selectionHN (vector<vector<vector<bool> > > X, vector<double> F){

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

    //vector<vector<vector<bool> > > Mum(P,vector<vector<bool> > (N,vector<bool>(2,false)));
    vector<vector<vector<bool> > > Mum(P);
    for(int i=0;i<P;i++){
        Mum[i].resize(N);
        for(int j=0;j<N;j++){
            Mum[i][j].resize(2,false);
        }
    }

    vector<vector<vector<bool> > > Dad(P);
    for(int i=0;i<P;i++){
        Dad[i].resize(N);
        for(int j=0;j<N;j++){
            Dad[i][j].resize(2,false);
        }
    }


    for(int p=0;p<P;p++){
        double r = rng->get()*cumSum;
        for(int q=0;q<P;q++){
            if ((lower[q]<=r) && (r<upper[q])){
                Mum[p] = X[q];
                break;
            }
        }
    }

    //vector<vector<vector<bool> > > Dad(P,vector<vector<bool> > (N,vector<bool>(2,false)));
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
        int crossOverGenetic = rng->get()*(N-1)+1;
        for(int i=crossOverGenetic;i<N;i++){
            X[p][i][0] = Dad[p][i][0];
        }

        int crossOverAdaptive = rng->get()*(N-1)+1;
        for(int i=crossOverAdaptive;i<N;i++){
            X[p][i][1] = Dad[p][i][1];
        }
    }
    return X;
}


int main (int argc, char **argv){

    if (argc < 3) { std::cerr << "\nUsage: ./hintonBatchVolatileConfigurable configfile logdir genomeMutationRate targetMutationRate numMutable seed(optional)"; return -1; }

    srand(time(NULL)); // note may not be different for simultaneously launched progs
    int seed = rand();
    if(argc==7){
        seed = std::stoi(argv[6]);
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

    int fitnessMode = conf.getInt("fitnessMode", 0);

    int G = conf.getInt("generations", 1000);
    int P = conf.getInt("populationSize", 1000);
    int T = conf.getInt("lifetimeLearningTrials", 1000); // NOTE SET T=-1 IN CONFIG TO USE PARAM LINE INPUT ARGV[4] TO SPECIFY IT INSTEAD
    int N = conf.getInt("numAlleles", 20);
    double initialQs = conf.getDouble("initialProportionLearnable", 0.5);

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

    double norm = 1.0/((double)P*(double)N);

    vector<int> rp= randPermute(N);
    vector<int> mutableID(numMutable);
    for(int i=0;i<numMutable;i++){
        mutableID[i]=rp[i];
    }

    double f_min = conf.getDouble("f_min", 1.0/(double)P);

    vector<bool> target(N,1);

    //
    //vector<vector<vector<bool> > > X(P,vector<vector<bool> > (N,vector<bool>(2,false)));

    vector<vector<vector<bool> > > X(P);
    for(int i=0;i<P;i++){
        X[i].resize(N);
        for(int j=0;j<N;j++){
            X[i][j].resize(2,false);
        }
    }


    for(int p=0;p<P;p++){
        for(int i=0;i<N;i++){
            X[p][i][0]=rng->get()<0.5;          // genetic state
            X[p][i][1]=rng->get()<initialQs;    // whether plastic
        }
    }

    vector<double> p0(G,0.0);
    vector<double> p1(G,0.0);
    vector<double> p2(G,0.0);

    // MAIN EVOLUTINOARY LOOP

    for(int g=0;g<G;g++){

        // extract states
        int c0=0;
        int c1=0;
        int c2=0;
        for(int p=0;p<P;p++){
            for(int i=0;i<N;i++){
                if(X[p][i][1]){
                    c2++;
                } else {
                    if(X[p][i][0]==target[i]){
                        c1++;
                    } else {
                        c0++;
                    }
                }
            }
        }
        p0[g]=(double)c0*norm;
        p1[g]=(double)c1*norm;
        p2[g]=(double)c2*norm;


        // vary target
         for(int i=0;i<numMutable;i++){
                if(rng->get()<targetMutationRate){
                    target[mutableID[i]] = !target[mutableID[i]];
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

                vector<vector<bool> > x(N/chunks,vector<bool> (2,false));
                vector<bool> t(N/chunks,false);
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
                    F[p] = f/(double)chunks;
                }
            } break;
        }

        X = selectionHN(X,F);


        // mutation
        if(mutationRate>0){
            for(int p=0;p<P;p++){
                for(int i=0;i<N;i++){
                    if(rng->get()<mutationRate){
                        X[p][i][0] = !X[p][i][0];
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

    return 0;
}



