//
// Class for sequences of times for genetic algorithm
// optimisation of period detection etc
//

#include <stdlib.h>
#include <iostream.h>
#include "trm/subs.h"

class sequence{
public:
  sequence() : nseq(0) {};
  sequence(int nsq);
  sequence(int nsq, int nmax, long int &seed);
  ~sequence();

  sequence(const sequence &obj);
  sequence &operator=(const sequence &head);
  int &operator[](int i);

  bool differ();
  void rmutate(int nmax, long int &seed);
  void lmutate(int nmax, int mdiff, long int &seed);
  void order();
  void ran(int nmax, long int &seed);
  void set_nseq(int nsq);

  friend std::ostream &operator<<(std::ostream &stream, const sequence &obj);  
  friend sequence splice(sequence &par1, sequence &par2, long int &seed);

private:
  int nseq;
  int *seq;
};

// constructor of the right size

sequence::sequence(int nsq) : nseq(nsq) {
  seq = new int[nseq];
}

// constructor of a non-repetitive random sequence

sequence::sequence(int nsq, int nmax, long int &seed) : nseq(nsq) {
  seq = new int[nseq];
  this->ran(nmax, seed);
}

// sets number in a sequence

void sequence::set_nseq(int nsq){
  if(nseq) delete seq;
  nseq  = nsq;
  seq = new int[nseq];
  return;
}

// destructor

sequence::~sequence(){
  if(nseq) delete seq;
}

// copy constructor

sequence::sequence(const sequence &obj){
  nseq = obj.nseq;
  seq  = new int[nseq];
  for(int i=0; i<nseq; i++){
    seq[i] = obj.seq[i];
  }
}

// assignment

sequence &sequence::operator=(const sequence &obj){
  if(obj.nseq){
    if(nseq != obj.nseq){
      if(nseq) delete seq;
      nseq = obj.nseq;
      seq  = new int[nseq];
    }
    for(int i=0; i<nseq; i++){
      seq[i] = obj.seq[i];
    }
    return *this;
  }else{
    std::cerr << "Error in sequence &sequence::operator=(sequence &obj)" << std::endl;
    exit(EXIT_FAILURE);
  }
}

// array

int &sequence::operator[](int i){
  if(nseq){
    return seq[i];
  }else{
    std::cerr << "Error in int &sequence::operator[](int i)" << std::endl;
    exit(EXIT_FAILURE);
  }
}

// sets a random sequence with no repetition

void sequence::ran(int nmax, long int &seed){
  if(nseq){
    bool ok;
    int trial, is;
    for(int ns=0; ns < nseq; ns++){
      ok = false;
      while(!ok){
	trial = int(nmax*ran2(seed));
	ok = true;
	for(is=0; is < nseq-1; is++){
	  if(trial == seq[is]){
	    ok = false;
	    break;
	  }
	}
      }
      seq[ns] = trial;
    }
  }else{
    std::cerr << "void sequence::ran(long int &seed) error" << std::endl;
    exit(EXIT_FAILURE);
  }
}

// order orders a sequence into ascending order

void sequence::order(){
  if(nseq){
    heapsort(seq,nseq);
  }else{
    std::cerr << "void sequence::order() error" << std::endl;
    exit(EXIT_FAILURE);
  }
}

// checks that all are different

bool sequence::differ(){
  if(nseq){
    int is, js;
    for(is=1; is < nseq; is++){
      for(js=0; js < is; js++){
	if(seq[js] == seq[is]) return false;
      }
    }
    return true;
  }else{
    std::cerr << "bool sequence::differ() error" << std::endl;
    exit(EXIT_FAILURE);
  }
}

// random mutation -- changes one completely at random

void sequence::rmutate(int nmax, long int &seed){
  if(nseq){
    int n, i, j = (int)(nseq*ran2(seed));
    bool ok = false;
    while(!ok){
      n = (int)(nmax*ran2(seed));
      ok = true;
      for( i=0; i<nseq; i++){
	if(n == seq[i]){
	  ok = false;
	  break;
	}
      }
    }
    seq[j] = n;
  }else{
    std::cerr << "void sequence::mutate(int nmax, long int &seed) error" << std::endl;
    exit(EXIT_FAILURE);
  }
}

// local mutation -- changes one to a nearby slot

void sequence::lmutate(int nmax, int mdiff, long int &seed){
  if(nseq){
    int n, i, j = (int)(nseq*ran2(seed));
    bool ok = false;
    while(!ok){
      n = -1;
      while(n < 0 || n >= nmax){
	n = seq[j]+(int)(mdiff*(2*ran2(seed)-1));
      }
      ok = true;
      for( i=0; i<nseq; i++){
	if(n == seq[i]){
	  ok = false;
	  break;
	}
      }
    }
    seq[j] = n;
  }else{
    std::cerr << "void sequence::mutate(int nmax, long int &seed) error" << std::endl;
    exit(EXIT_FAILURE);
  }
}

// friend function splice to generate a new sequence fromtwo parents.

sequence splice(sequence &par1, sequence &par2, long int &seed){
  if(par1.nseq > 1 && par1.nseq == par2.nseq){
    bool used1 = false, used2 = false;
    sequence temp(par1.nseq);
    /*
    int i, j = (int)((par1.nseq-2)*ran2(seed))+1;
    for(i=0;i<j;i++){
      temp.seq[i] = par1.seq[i];
    }
    for(i=j;i<par2.nseq;i++){
      temp.seq[i] = par2.seq[i];
    }
    */

    // select from each parent at random making sure
    // that each is used at least once

    for(int i=0; i < par1.nseq-1; i++){
      if(ran2(seed) < 0.5){
	temp.seq[i] = par1.seq[i];
	used1 = true;
      }else{
	temp.seq[i] = par2.seq[i];
	used2 = true;
      }
    }
    if(used1 && used2){
      if(ran2(seed) < 0.5){
	temp.seq[par1.nseq-1] = par1.seq[par1.nseq-1];
      }else{
	temp.seq[par1.nseq-1] = par2.seq[par1.nseq-1];
      }
    }else if(used1){
      temp.seq[par1.nseq-1] = par2.seq[par1.nseq-1];
    }else{
      temp.seq[par1.nseq-1] = par1.seq[par1.nseq-1];
    }
    return temp;
  }else{
    std::cerr << "sequence &splice(sequence par1, sequence par2) error" << std::endl;
    exit(EXIT_FAILURE);
  }
}

std::ostream &operator<<(std::ostream &stream, const sequence &obj){
  if(obj.nseq){
    stream << obj.seq[0];
    for(int i = 1; i < obj.nseq; i++){
      stream << " " << obj.seq[i];
    }
  }
  return stream;
}
