#include <iostream>
#include <string>
#include <assert.h>

//
// these are givens - let's assume this code works.
//

typedef unsigned long long HashType;

#define uniqify_rc(f, r) ((f) < (r) ? (f) : (r))

#define twobit_repr(ch) ((ch) == 'A' ? 0LL : \
                         (ch) == 'T' ? 1LL : \
                         (ch) == 'C' ? 2LL : 3LL)

#define twobit_comp(ch) ((ch) == 'A' ? 1LL : \
                         (ch) == 'T' ? 0LL : \
                         (ch) == 'C' ? 3LL : 2LL)

HashType _hash(const char * kmer, const char k, HashType& _h, HashType& _r)
{
  HashType h = 0, r = 0;

  h |= twobit_repr(kmer[0]);
  r |= twobit_comp(kmer[k-1]);

  for (int i = 1, j = k - 2; i < k; i++, j--) {
    h = h << 2;
    r = r << 2;

    h |= twobit_repr(kmer[i]);
    r |= twobit_comp(kmer[j]);
  }

  _h = h;
  _r = r;

  return uniqify_rc(h, r);
}

void increment_count(HashType kmer)
{
  // stub code, for demo purposes
  std::cout << kmer << "\n";
}

//
// this is the function we're going to be looking at today.
// it's an optimized version of a function that extracts all
// DNA substrings of length k, calculates their forward and
// reverse-complement hashes, chooses the lower of the two, and
// passes that value to increment_count.
//

// https://tinyurl.com/titusCode

// this function counts k-mers

class KmerIterator {
public:
    std::string input_string;
    const unsigned int k;
    unsigned int pos;
    HashType bitmask;
    HashType f, r;
    bool _initialized;

    KmerIterator(std::string inp, unsigned int ksize) :
        input_string(inp), k(ksize)
    {
        // initialize a 2*k-sized bitmask to all ones.
        for (unsigned char i = 0; i < k; i++) {
            bitmask = (bitmask << 2) | 3;
        }
        f = 0;
        r = 0;
        _initialized = false;
    }

    bool finished() {
        return (pos >= input_string.length());
    }

    HashType next()
    {
        HashType h;
        if (!_initialized) {
            h = _hash(input_string.c_str(), k, f, r);
            pos = k;
            _initialized = true;
            return h;
        }

        unsigned char ch = *(input_string.c_str() + pos);

        // bitshift left by 2, giving room for new 2bit hash
        f = f << 2;
        f |= twobit_repr(ch);
        f &= bitmask;

        // do the same for the reverse complement
        r = r >> 2;
        r |= (twobit_comp(ch) << (k*2 - 2));
        // choose a single canonical representation
        h = uniqify_rc(f, r);

        pos += 1;
        return h;
    }
};

void count_kmers(const std::string &input_string, const unsigned int k)
{
    KmerIterator ki(input_string, k);
    HashType h;

    while (!ki.finished()) {
        h = ki.next();
        increment_count(h);
    }
}

void slow_count_kmers(const std::string &s, const unsigned int k)
{
    HashType f, r;
    for (unsigned int i = 0; i < s.length() - k + 1; i++){
        HashType h = _hash(s.c_str() + i, k, f, r);
        increment_count(h);
    }
}


int main()
{
  count_kmers("ATGGGACCAGATAGAGCCAGAGGACACATTAGGACGAT", 3);
}
