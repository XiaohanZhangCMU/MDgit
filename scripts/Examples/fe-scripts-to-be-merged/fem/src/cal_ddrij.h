/*
 * cal_ddrij.h
 * by Keonwook Kang kwkang@lanl.gov 
 * Last Modified : Wed May 18 22:35:50 MDT 2011
 */

#define MAXBONDS 12000

//namespace MDPP_NS {

class ComputeDisregistry: public MDFrame
{
public:
    ComputeDisregistry();
    ~ComputeDisregistry();

    int nb;                  // Number of bonds
    double slip_normal[10];
    struct BONDS 
    { 
        int i, j, type;
        Vector3 center;
    } *bonds;
    //} bonds [MAXBONDS];
    struct DDRIJ
    {
        int tag;
        double norm;
        Vector3 disreg;
    } *ddrij_sav;

    void alloc();
    void initparser();
    void initvars();
    void calddrij();
    void plot_ddrij();
    int identify_bonds();
    void write_bonds(char *);
    void read_bonds(char *);
    void plot_bonds();
    void retag_ddrij();
    void write_ddrij(char *);
    void read_ddrij(char *);
    void potential();

    virtual int exec(const char *name);
    
    //private:

};

//}
