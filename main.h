/*
 *
 */

#ifndef MAIN_H_
#define MAIN_H_

#define PROGRAM "bamDCS"
#define VERSION "v2.4"
#define AUTHORS "yerui"
#define CONTACT "yerui@connect.hku.hk"
#define REMARKS "(generate double strand consensus reads for low-frequency mutation)"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <utility>
#include <algorithm>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "gzstream.h"
#include "compass_search.hpp"
#include "boost/math/distributions/chi_squared.hpp"

using namespace __gnu_cxx;
using namespace std;
using namespace SeqLib;

typedef unsigned long ulong;

class Option {
    public:
        Option():
            baseQuaCutoff(20),
            mapQuaCutoff(30),
            minSupOnEachStrand(3),
            maxSupOnEachStrand(3000),
            Ncutoff(0.1),
            minFractionInFam(0.002),
            freqPrecision(0.00001),
            lhrGapCutoff(5.0),
            phredOffset(33),
            minSupOnHaplotype(3),
            filtSoftClip(false),
            outBamFile(""),
            debug(false),
            pvalue(0.001),
            pcrError(1.0e-5),
            softEndTrim(5) {}

        ~Option(){}

        int baseQuaCutoff;
        int mapQuaCutoff;
        ulong minSupOnEachStrand;
        ulong maxSupOnEachStrand;

        double Ncutoff;
        double minFractionInFam;
        double freqPrecision;
        double lhrGapCutoff;

        int phredOffset;
        ulong minSupOnHaplotype;
        bool filtSoftClip;
        string outBamFile;
        bool debug;
        double pvalue;
        double pcrError;
        ulong softEndTrim;
};

typedef map<string, vector<SeqLib::BamRecord> > mStrBrV;
typedef map<char,double> mCharDouble;
typedef map<char, ulong> mCharUlong;
typedef map<string,ulong> mStrUlong;

typedef pair<char,ulong> pCharUlong;
typedef pair<double, set<char> > pDoubleCharSet;

typedef vector< set<char> > vCharSet;
typedef vector<string> vString;

// descending sort pair
inline bool _cmpByFirst(const pDoubleCharSet &a, const pDoubleCharSet &b) { return a.first > b.first; }
inline bool _cmpBySecond(const pCharUlong &a, const pCharUlong &b) { return a.second > b.second; }

inline double errorRate(const char &q, const Option &opt);
inline vector<double> quaToErrorRate(const string &qs, const Option &opt);
inline string _itoa( size_t &i );
inline ulong countN(const string &s);

inline bool lowQuality(const char &, const Option &);
inline string reverseComplement(const string &);

void usage();
SeqLib::BamHeader removeReadGroup( const SeqLib::BamHeader & );
void calibrateFam(mStrBrV &);
bool testCigarFam(mStrBrV &, mStrBrV &, const string &, const Option &);
void printConsensusRead(ogzstream &,ogzstream &,mStrBrV &,mStrBrV &,const string &,const Option &,const string &, SeqLib::BamWriter &);
string getQuaFromPvalue( const vector< mCharDouble > &quaV, const string &s, const Option &opt );
string adjust_p(const string &qs, const Option &opt);

//vCharSet hetPoint(const BamRecordVector &, const Option &);
vector< mCharDouble > hetPoint(const BamRecordVector &, const Option &);

//vCharSet zipHetPoint(const vCharSet &, const vCharSet &);
vector< mCharDouble > zipHetPoint(const vector< mCharDouble > &, const vector< mCharDouble > &, const Option &opt);
vString consensusSeq(const BamRecordVector &, const BamRecordVector &, const vector< mCharDouble > &, const Option &);

//set<char> llh_genotype(const string &, const string &, const Option &);
mCharDouble llh_genotype(const string &, const string &, const Option &);

pDoubleCharSet maxLogLikelihood(const string &, const vector<double> &, const vector<pCharUlong> &, const Option &, const int);
double minus_llh_3nt( int m, double x[], const vector<pCharUlong> &v, const string &s, const vector<double> &e,  double lower[], double upper[], double sumBound );
double calculate_llh(const string &, const vector<double> &, mCharDouble &);
string ignoreError(const string &, const vector< mCharDouble > &);
void trimEnd(SeqLib::BamRecord &br, const Option &opt);
//SeqLib::Cigar trimCigar(SeqLib::BamRecord &br, const SeqLib::Cigar &cg, const Option &opt);
size_t getGenomePosition(size_t pos, ulong i, const Cigar &cg);

// sub functions

void usage() {
    cout << "Program: " << PROGRAM << "  " << REMARKS << "\n"
        "Version: " << VERSION << "\n"
        "Authors: " << AUTHORS << "\n"
        "Contact: " << CONTACT << "\n\n"

        "Options: " << PROGRAM << " in.bam out_prefix\n\n"

        "    -q [i]     base quality cutoff [20]\n"
        "    -Q [i]     map quality cutoff [30]\n"

        "    -s [i]     min support num on each strand [3]\n"
        "    -S [i]     max support num on each strand [3000]\n"

        "    -N [f]     max fraction of N on consensus read [0.1]\n"
        "    -f [f]     min fraction of alterative allele in a read family [0.002]\n"

        "    -e [f]     precision of allele frequency while calculate likelihood ratio [0.00001]\n"
        "    -g [f]     gap between likelihood ratios of major and rest genotypes [5.0]\n"
        "    -x [i]     Encoding offset for phred quality scores [33]\n"
        "    -t [i]     min support num to construct a haplotype seq [3]\n"
        "    -c         discard soft-clipping reads [false]\n"
        "    -C [i]     soft trim N base from both ends of read [5]\n"

        "    -o [s]     output bam File directly []\n"
        "    -d         debug mode [false]\n"
        "    -h         help\n"
        "    -v         version\n"
        "\n";
}

SeqLib::BamHeader removeReadGroup( const SeqLib::BamHeader &oldHead )
{
    istringstream itm( oldHead.AsString() );
    string line;
    string s("");

    while (getline(itm, line) ) {
        size_t i = line.find("@RG");

        if ( i == 0 ) {
            s += "@RG\tID:foo\tSM:bar\n";
        }
        if ( i != 0 ) {
            s += line + "\n";
        }
    }

    // add rg if the old head does not have rg
    if ( s.find("@RG") == string::npos ) {
        s += "@RG\tID:foo\tSM:bar\n";
    }

    SeqLib::BamHeader newHead(s);
    return newHead;
}

inline double errorRate(const char &q, const Option &opt)
{
    return pow( 10, (double)(opt.phredOffset-(int)q)/10 );
}

inline vector<double> quaToErrorRate(const string &qs, const Option &opt)
{
    vector<double> eV;
    for ( size_t i(0); i != qs.size(); i++ )  eV.push_back( errorRate(qs[i], opt) );
    return eV;
}

inline string _itoa( size_t &i )
{
    ostringstream osm;
    osm << i;
    return osm.str();
}

inline ulong countN(const string &s)
{
    ulong n(0);
    for ( size_t i(0); i != s.size(); i++ ) {
        if ( s[i] == 'N' ) n++;
    }
    return n;
}

inline bool lowQuality(const char &q, const Option &opt)
{
    return ( int(q) - opt.phredOffset < opt.baseQuaCutoff );
}

inline string reverseComplement(const string &seq)
{
    string str("");
    for ( string::const_reverse_iterator it = seq.rbegin(); it != seq.rend(); it++) {
        switch (*it) {
            case 'A': str += "T";    break;
            case 'C': str += "G";    break;
            case 'G': str += "C";    break;
            case 'T': str += "A";    break;
            case 'N':
            default:  str += "N";    break;
        }
    }

    return str;
}

inline char errorRateToChar(const double &q, const Option &opt)
{
    return char( opt.phredOffset - (int)(-10 * log(q) / log(10.0) ) );
}

void calibrateFam(mStrBrV &)
{
    return;
}

// check family size on cigar
bool testCigarFam(mStrBrV &watsonFam, mStrBrV &crickFam, const string &cg, const Option &opt)
{
    mStrBrV::const_iterator w = watsonFam.find( cg );
    mStrBrV::const_iterator c =  crickFam.find( cg );

    if (    w != watsonFam.end()
            && c !=  crickFam.end()
            && w->second.size() >= opt.minSupOnEachStrand * 2
            && c->second.size() >= opt.minSupOnEachStrand * 2
            && w->second.size() <= opt.maxSupOnEachStrand * 2
            && c->second.size() <= opt.maxSupOnEachStrand * 2
       )
        return true;
    else
        return false;
}

// output consensus Read pair into fq files
void printConsensusRead(
        ogzstream & fq1,
        ogzstream & fq2,
        mStrBrV & watsonFam,
        mStrBrV & crickFam,
        const string & cg,
        const Option & opt,
        const string & chrBegEnd,
        SeqLib::BamWriter & writer
        ) {
    // length of read1 & read2 are same, So connect them to simplify workflow

    // 0-based index of heterozygous point --> vector of allele set
//    vCharSet wHetPos, cHetPos, sameHetPos;
    vector< mCharDouble > wHetPos, cHetPos, sameHetPos;

    // find heterzygous point on each fam
    wHetPos = hetPoint(watsonFam[cg], opt);
    cHetPos = hetPoint( crickFam[cg], opt);

    // find consistent het point by comparing watson & crick family
    sameHetPos = zipHetPoint(wHetPos, cHetPos, opt);

    // watson & crick should be concordant on hom point. if not, set N
    vString seqV = consensusSeq(watsonFam[cg], crickFam[cg], sameHetPos, opt);

    string Qname("@");
    Qname += chrBegEnd + ":" + cg + ":";

    if ( seqV.empty() ) return;

    for ( size_t i(0); i != seqV.size(); i++ ) {
        string seq = seqV[i];
        size_t length = seq.size()/2;
        string quaStr = getQuaFromPvalue( sameHetPos, seq, opt );

        string id1 = Qname + _itoa(i) + "/1";
        string id2 = Qname + _itoa(i) + "/2";

        string rd1 = seq.substr( 0, length );
        string rd2 = reverseComplement( seq.substr(length) );

        string quaStr1 = quaStr.substr( 0, length );
        string quaStr2 = quaStr.substr( length    );

        reverse( quaStr2.begin(), quaStr2.end() );

        /*
        for ( size_t j(0); j != length; j++ ) {
            quaStr1 += ( rd1[j] == 'N' ? '$' : 'J' );
            quaStr2 += ( rd2[j] == 'N' ? '$' : 'J' );
        }
*/
//   fq1 << id1 << "\n" << rd1 << "\n+\n" << quaStr1 << endl;
  //      fq2 << id2 << "\n" << rd2 << "\n+\n" << quaStr2 << endl;

        if ( opt.outBamFile.size() > 0 ) {
            SeqLib::BamRecord br1 = watsonFam[cg].at(0);
            SeqLib::BamRecord br2 = watsonFam[cg].at(1);

//            br1.SetQname( chrBegEnd + ":" + cg );
//            br2.SetQname( chrBegEnd + ":" + cg );
            if ( opt.softEndTrim > 0 ) {
     //           br1.SetCigar( trimCigar(br1, br1.GetCigar(), opt) );
       //         br2.SetCigar( trimCigar(br2, br2.GetCigar(), opt) );
                trimEnd(br1, opt);
                trimEnd(br2, opt);
            }

            br1.SetSequence( rd1 );
            br2.SetSequence( seq.substr(length) );

            reverse( quaStr2.begin(), quaStr2.end() );
            br1.SetQualities( quaStr1, 33 );
            br2.SetQualities( quaStr2, 33 );

            br1.RemoveAllTags();
            br2.RemoveAllTags();

            br1.AddZTag("RG", "foo");
            br2.AddZTag("RG", "foo");

            writer.WriteRecord( br1 );
            writer.WriteRecord( br2 );
        }
    }
}

// get heterozygous points based on watson or crick family only
vector< mCharDouble > hetPoint(const BamRecordVector &brV, const Option &opt)
{
//    vCharSet pt;
    vector< mCharDouble > pt;
    vString seqV, quaV;

    // too few reads to support heterozygous point. two reads form a pair.
    if ( brV.size() < opt.minSupOnEachStrand * 2 )
        return pt;

    string seq(""), qua("");
    for ( auto & br : brV ) {
        if ( !br.ReverseFlag() ) { // + strand
            seq = br.Sequence();
            qua = br.Qualities();
        }
        else { // - strand
            seq += br.Sequence();
            qua += br.Qualities();

            seqV.push_back(seq);
            quaV.push_back(qua);
        }
    }

    size_t length = seq.size();

    // fetch each column of alleles and quality scores
    for ( size_t j(0); j != length; j++ ) {
        string base(""), qual("");

        for ( size_t i(0); i != seqV.size(); i++ ) {
            base += seqV[i][j];
            qual += quaV[i][j];
        }

//        string adj_qua = adjust_p(qual, opt);
        pt.push_back( llh_genotype(base, qual, opt) );
    }

    return pt;
}

string getQuaFromPvalue( const vector< mCharDouble > &quaV, const string &s, const Option &opt )
{
    ostringstream q("");

    for ( size_t i(0); i != s.size(); i++ ) {

        if ( s[i] == 'N' ) {
            q << '$';
        }
        else {
            mCharDouble::const_iterator it = quaV[i].find( s[i] );

            if ( it != quaV[i].end() ) {
                double f = -10.0 * log( it->second )/log(10.0);

       //         if ( f > 41.0 ) f = 41.0;
         //       if ( f < 0.0  ) f = 0.0;
                q << (char)( opt.phredOffset + int(f) );
            }
            else {
                q << '$';
            }

        }
    }

    return q.str();
}

vector< mCharDouble > zipHetPoint(const vector< mCharDouble > &w, const vector< mCharDouble > &c, const Option &opt)
{
//    vCharSet samePt;
    vector< mCharDouble > samePt;

    if ( w.empty() || c.empty() )
        return samePt;

    for ( size_t j(0); j != w.size(); j++ ) {
        //set<char> nt;
        mCharDouble nt;

        for ( mCharDouble::const_iterator wi = w[j].begin(); wi != w[j].end(); wi++ ) {
            if ( wi->second > opt.pvalue ) continue;

            mCharDouble::const_iterator ci = c[j].find( wi->first );
            if ( ci != c[j].end() ) {
                if ( ci->second > opt.pvalue ) continue;
               // nt.insert( *it );

                nt[ wi->first ] = wi->second + ci->second - wi->second * ci->second
                    + 10 * pow(opt.pcrError,2);
            }
        }

        samePt.push_back(nt);
    }

    return samePt;
}

//set<char> llh_genotype(const string &s, const string &q, const Option &opt)
mCharDouble llh_genotype(const string &s, const string &q, const Option &opt)
{
    // allele set which will be returned
//    set<char> ntS;
    mCharDouble ntP; // nt --> pvalue
    boost::math::chi_squared X2_dist(1);

    // check frequent of alleles
    mCharUlong fr;

    double depth( s.size() );
    double small_diff(1e-10);

    for ( size_t i(0); i != depth; i++ ) {
        if ( lowQuality(q[i], opt) || s[i] == 'N' )  continue;
//        if ( s[i] == 'N' )  continue;
        fr[ s[i] ]++;
    }

    if ( fr.empty() ) {
//        return ntS;
        return ntP;
    }

    // sort by frequent
    vector<pCharUlong> ntV( fr.begin(), fr.end() );
    vector<double> errV = quaToErrorRate(q, opt);

    // only one allele
    if ( ntV.size() == 1 ) {
       // if ( ntV[0].second >= opt.minSupOnEachStrand ) {
         //   ntS.insert( ntV[0].first );
        //}

        pDoubleCharSet tmp = maxLogLikelihood(s,errV,ntV,opt,1);
        for ( auto &e : errV ) tmp.first -= ( log(e/3) ); // null hypothesis

        if ( tmp.first <= 0 ) {
     //       cout << "only one: " << tmp.first << ' ' << ntV[0].first << ":" << ntV[0].second << endl;
            tmp.first += small_diff;
        }

        ntP[ ntV[0].first ] = 1 - boost::math::cdf(X2_dist, 2*tmp.first);
        //return ntS;
        return ntP;
    }

//    if ( ntV.size() > 1 )
  //      sort( ntV.begin(), ntV.end(), _cmpBySecond ); // descending sort


    if ( ntV.size() == 2 ) {
        pDoubleCharSet two = maxLogLikelihood(s,errV,ntV,opt,2);

        vector<pCharUlong> ntV_1, ntV_2;
        ntV_1.push_back( ntV[1] );
        ntV_2.push_back( ntV[0] );

        pDoubleCharSet t1 = maxLogLikelihood(s,errV,ntV_1,opt,1);
        pDoubleCharSet t2 = maxLogLikelihood(s,errV,ntV_2,opt,1);

        if ( two.first - t1.first <= 0 ) {
//            cout << "only two t1: " << t1.first << ' ' << two.first << endl;
            t1.first -= small_diff;
        }

        if ( two.first - t2.first <= 0 ) {
  //          cout << "only two t2: " << t2.first << ' ' << two.first << endl;
            t2.first -= small_diff;
        }


        ntP[ ntV[0].first ] = 1 - boost::math::cdf(X2_dist, 2*(two.first - t1.first) );
        ntP[ ntV[1].first ] = 1 - boost::math::cdf(X2_dist, 2*(two.first - t2.first) );

        return ntP;
    }

    sort( ntV.begin(), ntV.end(), _cmpBySecond ); // descending sort

    if ( ntV.size() > 3 )
        ntV.pop_back();

    if ( opt.debug )
        for ( auto &nt : ntV ) cout << nt.first << "=>" << nt.second << ' ';

    pDoubleCharSet three = maxLogLikelihood(s,errV,ntV,opt,3);

    vector<pCharUlong> ntV1, ntV2, ntV3;
    ntV1.push_back( ntV[1] );
    ntV1.push_back( ntV[2] );

    ntV2.push_back( ntV[0] );
    ntV2.push_back( ntV[2] );

    ntV3.push_back( ntV[0] );
    ntV3.push_back( ntV[1] );

    pDoubleCharSet tm1 = maxLogLikelihood(s,errV,ntV1,opt,2);
    pDoubleCharSet tm2 = maxLogLikelihood(s,errV,ntV2,opt,2);
    pDoubleCharSet tm3 = maxLogLikelihood(s,errV,ntV3,opt,2);

    if ( three.first - tm1.first <= 0 ) {
        //cout << "three tm1: " << tm1.first << endl;
        tm1.first -= small_diff;
    }

    if ( three.first - tm2.first <= 0 ) {
        //cout << "three tm2: " << tm2.first << endl;
        tm2.first -= small_diff;
    }

    if ( three.first - tm3.first <= 0 ) {
        //cout << "three tm3: " << tm3.first << endl;
        tm3.first -= small_diff;
    }

    ntP[ ntV[0].first ] = 1 - boost::math::cdf(X2_dist, 2*(three.first - tm1.first) );
    ntP[ ntV[1].first ] = 1 - boost::math::cdf(X2_dist, 2*(three.first - tm2.first) );
    ntP[ ntV[2].first ] = 1 - boost::math::cdf(X2_dist, 2*(three.first - tm3.first) );

    return ntP;





    /*
    // delete pair<allele, supportNum> which has too few support reads or too small fraction
    while ( ntV.size() ) {
        if ( ntV.back().second < opt.minSupOnEachStrand || ntV.back().second / depth < opt.minFractionInFam )
            ntV.pop_back();
        else
            break;
    }

    // none allele remain
    if ( ntV.empty() ) {
  //      return ntS;
        return ntP;
    }

    // only one allele
    if ( ntV.size() == 1 ) {
        if ( ntV[0].second >= opt.minSupOnEachStrand )
            ntS.insert( ntV[0].first );

        return ntS;
    }

    // two or more alleles
    vector<double> errV = quaToErrorRate(q, opt);
    vector<pDoubleCharSet> llhV;

    if ( ntV.size() > 1 ) {
        llhV.push_back( maxLogLikelihood(s,errV,ntV,opt,1) );
        llhV.push_back( maxLogLikelihood(s,errV,ntV,opt,2) );

        if ( ntV.size() > 2 ) {
            llhV.push_back(  maxLogLikelihood(s,errV,ntV,opt,3) );
        }
    }

    for ( auto & p : llhV ) {
        if ( p.second.size() == 3 ) p.first -= 2 * opt.lhrGapCutoff;
        if ( p.second.size() == 2 ) p.first -=     opt.lhrGapCutoff;
    }

    if ( llhV.size() > 1 )
        sort(llhV.begin(), llhV.end(), _cmpByFirst); // desending sort

    // four alleles is really rare, so do not consider it here
    if ( llhV.empty() ) {
        return ntS;
    }
    else
        return llhV[0].second;

        */
}

pDoubleCharSet maxLogLikelihood(const string &s, const vector<double> &e, const vector<pCharUlong> &v, const Option &opt, const int mode)
{
    double llh(0.0);
    set<char> cSet;

    if ( mode == 1 ) {
        char a( v[0].first );
        for ( size_t i(0); i != s.size(); i++ ) {
            if ( s[i] == a )
                llh += log( 1 - e[i] );
            else
                llh += log( e[i]/3 );
        }

        cSet.insert(a);
        return make_pair(llh, cSet);
    }
    else if ( mode == 2 ) {
        //initial value
        char a( v[0].first ), b( v[1].first );
        double total( v[0].second + v[1].second );
        double af( v[0].second / total);
        double bf( 1 - af );

        for ( size_t i(0); i != s.size(); i++ ) {
            if ( s[i] == 'N' ) continue;

                  s[i] == a ? ( llh += log( (1 - e[i] * 4/3) * af + e[i]/3 ) )
                : s[i] == b ? ( llh += log( (1 - e[i] * 4/3) * bf + e[i]/3 ) )
                :             ( llh += log( e[i] / 3 )                       )
                ;
        }

        double mllh = llh; // max log likelihood
        double oldaf = af;
        long step(1);

        if (opt.debug) cout << "~org af bf llh: " << af << ' ' << bf << ' ' << llh << endl;

        // calculate max likelihood step by step
        while (1) {
            af += opt.freqPrecision * step;
            if ( af > 1 ) af = 1;
            if ( af < 0 ) af = 0;
            bf = 1 - af;

            llh = 0;
            for ( size_t i(0); i != s.size(); i++ ) {
                if ( s[i] == 'N' ) continue;

                      s[i] == a ? ( llh += log( (1 - e[i] * 4/3) * af + e[i]/3 ) )
                    : s[i] == b ? ( llh += log( (1 - e[i] * 4/3) * bf + e[i]/3 ) )
                    :             ( llh += log( e[i]/3 )                         )
                    ;
            }

            if (opt.debug) cout << "~loop af bf step llh: " << af << ' ' << bf << ' ' << step << ' ' << llh << endl;

            if ( llh > mllh ) {
                mllh = llh;
                step *= 2;  // the direction is right, increase step length
                oldaf = af;
                if ( af == 1 ) step = -1;
                if ( af == 0 ) step = 1;
            }
            else {
                af = oldaf;
                if ( step == 1 )
                    step = -1;
                else if ( step == -1) {
                    bf = 1 - af;
                    break;
                }
                else
                    step = (step > 0 ? 1 : -1);
            }
        }

        if ( af != 0 )   cSet.insert(a);
        if ( bf != 0 )   cSet.insert(b);

        return make_pair(mllh, cSet);
    }
    else if ( mode == 3 ) {
//
        double total(0.0);
//        mCharDouble mBaseFreq;
        for ( size_t i(0); i != 3; i++ ) {
            total += v[i].second;
        }
//
//       for ( size_t i(0); i != 3; i++ )
//           mBaseFreq[ v[i].first ] = v[i].second/total;
//
//       llh = calculate_llh(s, e, mBaseFreq);

        // although 3 alleles, only two free variables, maximize the function of two variables.
        // future plan

        double delta(0.01);
        double delta_tol(1e-8);
        double fx;
        int k;
        int k_max(20000);
        int m = 2;

        double *x;
        double *x0;
        double *lower;
        double *upper;

        lower = new double[m];
        upper = new double[m];
        lower[0] = 0.0;
        upper[0] = 1.0;
        lower[1] = 0.0;
        upper[1] = 1.0;

        x0 = new double[m];
        if (opt.debug) cout << "  Test COMPASS_SEARCH with the function.\n";

        x0[0] = v[0].second/total;
        x0[1] = v[1].second/total;

        if (opt.debug) {
            r8vec_print ( m, x0, "  Initial point X0:" );
            cout << "\n";
            cout << "  F(X0) = " << minus_llh_3nt( m, x0, v, s, e, lower, upper, 1.0 ) << "\n";
        }

        x = compass_search ( minus_llh_3nt, m, x0, v, s, e, lower, upper, 1.0, delta_tol, delta, k_max, fx, k );

        if (opt.debug) {
            r8vec_print ( m, x, "  Estimated minimizer X1:" );
            cout << "\n";
            cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
        }

        // first two allele
        for ( int i(0); i != 2; i++ ) {
            if ( x[i] >= opt.freqPrecision  )   cSet.insert(v[i].first);
        }

        // the third allele
        if ( 1 - x[0] - x[1] >= opt.freqPrecision )   cSet.insert(v[2].first);

        delete [] lower;
        delete [] upper;
        delete [] x;
        delete [] x0;

        return make_pair(-fx,cSet); // -fx is the biggest llh
    }
    else {
        cerr << "only support mode 1,2,3" << endl, exit(1);
    }
}

double calculate_llh(const string &s, const vector<double> &e, mCharDouble &mBaseFreq)
{
    double llh(0);
    if ( s.empty() || e.empty() )
        cerr << "error occur when calculate_llh " << endl, exit(1);

    for ( size_t i(0); i != s.size(); i++ ) {
        if ( s[i] == 'N' ) continue;

        mCharDouble::iterator it = mBaseFreq.find( s[i] );
        if ( it != mBaseFreq.end() ) {
            llh += log( (1 - e[i] * 4/3) * (it->second) + e[i]/3 );
        }
        else {
            llh += log( e[i]/3 );
        }
    }

    return llh;
}

// minimizing -llh equals to maximize llh
double minus_llh_3nt( int m, double x[], const vector<pCharUlong> &v, const string &s, const vector<double> &e, double lower[], double upper[], double sumBound )
{
    double llh(0.0);

    for ( int i(0); i != m; i++ ) {
        if ( x[i] < lower[i] ) x[i] = lower[i];
        if ( x[i] > upper[i] ) x[i] = upper[i];
    }

    if ( x[0] + x[1] > sumBound ) x[1] = sumBound - x[0];

    for ( size_t i(0); i != s.size(); i++ ) {
        if ( s[i] == 'N' ) continue;

              s[i] == v[0].first ? ( llh += log( (1 - e[i] * 4/3) * x[0] + e[i]/3 ) )
            : s[i] == v[1].first ? ( llh += log( (1 - e[i] * 4/3) * x[1] + e[i]/3 ) )
            : s[i] == v[2].first ? ( llh += log( (1 - e[i] * 4/3) * ( 1-x[0]-x[1] ) + e[i]/3 ) )
            :                      ( llh += log( e[i] / 3 )                       )
            ;
    }

    return -llh;
}

vString consensusSeq(const BamRecordVector &w, const BamRecordVector &c, const vector< mCharDouble > &sameHetPos, const Option &opt)
{
    vString seqV;
    if ( sameHetPos.empty() ) return seqV;

    mStrUlong mSeqN_w, mSeqN_c;
    string seq("");

    // count N
    size_t Ncnt(0);
    for ( auto & p : sameHetPos ) {
        if ( p.empty() )
            Ncnt++;
    }

    if ( Ncnt > sameHetPos.size() * opt.Ncutoff )
        return seqV;

    // get potential seq based on original read and heterozygous site
    for ( auto &br : w ) {
        if ( !br.ReverseFlag() ) { // + strand
            seq = br.Sequence();
        }
        else { // - strand
            seq += br.Sequence();
            mSeqN_w[ ignoreError(seq, sameHetPos) ]++;
        }
    }

    for ( auto &br : c ) {
        if ( !br.ReverseFlag() ) { // + strand
            seq = br.Sequence();
        }
        else { // - strand
            seq += br.Sequence();
            mSeqN_c[ ignoreError(seq, sameHetPos) ]++;
        }
    }

    for ( mStrUlong::iterator it = mSeqN_w.begin(); it != mSeqN_w.end(); it++ ) {
        if ( it->second >= opt.minSupOnHaplotype ) {
            size_t len = it->first.size() / 2;
            if (   countN( it->first.substr(  0,len) ) > opt.Ncutoff * len
                    && countN( it->first.substr(len,len) ) > opt.Ncutoff * len
               )
                continue;

            mStrUlong::iterator ct = mSeqN_c.find( it->first );

            if (   ct != mSeqN_c.end() && ct->second >= opt.minSupOnHaplotype )
                seqV.push_back( it->first );
        }
    }

    return seqV;
}

string ignoreError(const string &s, const vector< mCharDouble > &v)
{
    string seq("");
    for ( size_t i(0); i != s.size(); i++ ) {

        if ( v[i].empty() ) {
            seq += "N";
        }
        else if ( v[i].size() == 1 ) {
            seq += v[i].begin()->first;
        }
        else {

            mCharDouble::const_iterator it = v[i].find( s[i] );
            if ( it != v[i].end() ) {
                seq += s[i];
            }
            else {
                seq += "N";
            }
        }
    }
    return seq;
}

string adjust_p(const string &qs, const Option &opt)
{
    ostringstream o;

    map<double, vector<char> > m;
    vector<double> v = quaToErrorRate(qs, opt);

    if ( v.size() > 1000 ) { 
        sort( v.begin(), v.end() );

        ulong total( v.size() );
        for ( size_t i(0); i != total; i++ ) {
            m[ v[i] ].push_back( errorRateToChar(v[i], opt) );
        }

        for ( auto & q : qs) {
            double e = errorRate(q, opt);
            o << m[ e ].back();
            m[ e ].pop_back();
        }

        return o.str();
    }
    else {
        return qs;
    }
}

void trimEnd(SeqLib::BamRecord &br, const Option &opt)
{
    CigarField sc('S', opt.softEndTrim );
    Cigar cg = br.GetCigar();

    if ( cg.size() == 0 ) {
        return;
    }

    size_t head(0), tail(0);
    string ty("=XMIS");

    //trim head
    Cigar nc;

    //cerr << "old: " << cg << ' ';

    for (size_t i(0); i != cg.size(); i++ ) {
        if ( ty.find( cg[i].Type() ) != string::npos ) {
            head += cg[i].Length();
        }

        //cerr << head << "h " << endl;

        if ( head >= opt.softEndTrim ) {
            nc.add(sc);
            //cerr << br.Position() <<  ' ';

            br.SetPosition( getGenomePosition(br.Position(), opt.softEndTrim, cg) );

            //cerr << br.Position() << endl;

            size_t remain = head - opt.softEndTrim;
            if ( remain > 0 ) {
                CigarField tmp(cg[i].Type(), remain);
                nc.add(tmp);

                for (size_t j(i+1); j != cg.size(); j++) nc.add(cg[j]);
            }
            else {
                for (size_t j(i+1); j != cg.size(); j++) nc.add(cg[j]);
            }
            //cerr << "newR: " << nc << ' ';

            break;
        }
    }

    //trim tail
    Cigar rc;

    for ( int i( nc.size()-1 ); i >= 0; i-- ) {
        if ( ty.find( nc[i].Type() ) != string::npos ) {
            tail += nc[i].Length();
        }

        if ( tail >= opt.softEndTrim ) {
            rc.add(sc);

            int remain = tail - opt.softEndTrim;
            if ( remain > 0 ) {
                CigarField tmp(nc[i].Type(), remain);
                rc.add(tmp);
                for (int j(i-1); j >= 0; j--) rc.add(nc[j]);
            }
            else {
                for (int j(i-1); j >= 0; j--) rc.add(nc[j]);
            }
            //cerr << " new: " << rc << ' ';

            break;
        }
    }

    // reverse
    Cigar rev;

    for ( int i( rc.size()-1 ); i >= 0; i-- ) {
        rev.add( rc[i] );
    }

    if ( rev[1].Type() == 'D' ) {
        br.SetPosition( br.Position() - rev[1].Length() );
    }

    //cerr << " final: " << rev << endl;

    br.SetCigar( rev );

}
/*
SeqLib::Cigar trimCigar(SeqLib::BamRecord &br, const SeqLib::Cigar &cg, const Option &opt)
{
    CigarField sc('S', opt.softEndTrim );
    Cigar nc;

    if ( cg.size() == 0 ) {
        return cg;
    }
    if ( cg.size() == 1 ) {
    //    cout << br.Position() << ' ' << opt.softEndTrim << ' ';

        if ( cg.front().Type() == 'M' ) {
            nc.add(sc);
            CigarField md('M', cg.front().Length() - 2 * opt.softEndTrim );
            br.SetPosition( br.Position() + opt.softEndTrim );
            nc.add(md);
            nc.add(sc);
        }
        else {
            cerr << "bad cigar: " << cg << endl;
            exit(1);
        }
    }
    else if ( cg.size() == 2 ) {
  //      cout << br.Position() << ' ' << opt.softEndTrim << ' ';

        if ( cg.front().Type() == 'S' && cg.back().Type() == 'M' ) {
            if ( cg.front().Length() >= opt.softEndTrim ) {
                nc.add( cg.front() );
                CigarField md('M', cg.back().Length() - opt.softEndTrim );
                nc.add(md);
                nc.add(sc);
            }
            else {
                nc.add( sc );
                CigarField md('M', cg.back().Length() - opt.softEndTrim*2 + cg.front().Length() );
                br.SetPosition( br.Position() + opt.softEndTrim - cg.front().Length() );
                nc.add(md);
                nc.add(sc);
            }
        }
        else if ( cg.front().Type() == 'M' && cg.back().Type() == 'S' ) {
            if ( cg.back().Length() >= opt.softEndTrim ) {
                nc.add(sc);
                CigarField md('M', cg.back().Length() - opt.softEndTrim );
                br.SetPosition( br.Position() + opt.softEndTrim );
                nc.add(md);
                nc.add( cg.back() );
            }
            else {
                nc.add(sc);
                CigarField md('M', cg.front().Length() - opt.softEndTrim*2 + cg.back().Length() );
                br.SetPosition( br.Position() + opt.softEndTrim );
                nc.add(md);
                nc.add(sc);
            }
            }
            else {
            return cg;
            }
            }
            else {
            return cg;
            }

//    cout << br.Position() << endl;

return nc;
}
 */
size_t getGenomePosition(size_t pos, ulong i, const Cigar &cg)
{
    for ( vector<CigarField>::const_iterator it = cg.begin(); it != cg.end(); it++ ) { 
        char t = it->Type();
        size_t n = it->Length();

        switch (t) {
            case '=':
            case 'X':
            case 'M':
                if (i < n)  return (pos+i);
                else        (pos+=n, i-=n);
                break;
            case 'I':
                if (i < n)  return pos;
                else        (i -= n); 
                break;
            case 'S':
                if (i < n)  return -1; 
                else        i -= n;
                break;

            case 'N':
            case 'D':       pos += n;  break;
            case 'H':
            case 'P':                  break;

            default:  cerr << "unknown cigar field: " << cg << endl, exit(1);
        }   
    }   

    return pos + i;
}

#endif // MAIN_H_
