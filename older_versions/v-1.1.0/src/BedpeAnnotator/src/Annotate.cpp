#include "Annotate.h"

int getRefExonSeq(Reference &ref, string chr, uint32_t aa, uint32_t bb, int strand, vector<char> & seq)
{

    uint32_t refaa=ref.to_ref_pos(chr,aa);
    uint32_t refbb=ref.to_ref_pos(chr,bb);
    if(strand==0)
    {
        for(uint32_t x=refaa;x<=refbb;x++)
            seq.push_back(ref.getRefChar(x));
    }
    else
    {
        for(uint32_t x=refbb;x>=refaa;x--)
            seq.push_back(getCharComp(ref.getRefChar(x)));
    }

    return 0;
}

// have repeats
vector<fusion_junction_t> fjtvec;

int Annotate::assignJunctions(Gene& g, Reference & ref, Bedpe & bdpe) {

    fjtvec.clear();

    int index=0;

    for(int i=0;i<bdpe.size();i++)
    {
        bedpe_t b=bdpe.getBedpe(i);

        uint32_t bk5,bk3;
        bdpe.getPos(b,bk5,bk3);
        
        vector<int> gids5p;
        g.isInGene(b.chr1,bk5,gids5p);
        vector<int> gids3p;
        g.isInGene(b.chr2,bk3,gids3p);
        vector<junction_t> j1;
        g.getBestExon3(gids5p, bk5, b.strand1.compare("+")==0?0:1, j1);

        vector<junction_t> j2;
        g.getBestExon3(gids3p, bk3, b.strand2.compare("+")==0?1:0, j2);
        for (int x=0;x<j1.size(); x++)
        {
            for(int y=0;y<j2.size();y++)
            {
                fusion_junction_t fjt;
                fjt.fusion_id=index+1;

                if(j1[x].is5p==1 && j2[y].is5p==0)
                {
                    getRefExonSeq(ref,j1[x].chr,j1[x].pos1+1,j1[x].pos2,j1[x].strand,fjt.seq1);
                    getRefExonSeq(ref,j2[y].chr,j2[y].pos1+1,j2[y].pos2,j2[y].strand,fjt.seq2);
              
                    fjt.p5=j1[x];
                    fjt.p3=j2[y];
                    // to real position not 0 based
                    if(b.strand1.compare("-")==0)
                        fjt.p5.pos1=fjt.p5.pos1+1;
                    if(b.strand2.compare("+")==0)
                        fjt.p3.pos1=fjt.p3.pos1+1;

                }
                else if(j1[x].is5p==0 && j2[y].is5p==1)
                {
                    getRefExonSeq(ref,j2[y].chr,j2[y].pos1+1,j2[y].pos2,j2[y].strand,fjt.seq1);
                    getRefExonSeq(ref,j1[x].chr,j1[x].pos1+1,j1[x].pos2,j1[x].strand,fjt.seq2);
              
                    fjt.p5=j2[y];
                    fjt.p3=j1[x];
                    // to real position not 0 based
                    if(b.strand2.compare("-")==0)
                        fjt.p5.pos1=fjt.p5.pos1+1;
                    if(b.strand1.compare("+")==0)
                        fjt.p3.pos1=fjt.p3.pos1+1;

                }
                fjtvec.push_back(fjt);

            }
        }
        index++;
    }

    return 0;
}


int getInframe(fusion_junction_t & ft)
{
    if (ft.p5.isCoding==1 && ft.p3.isCoding==1 && ft.p5.coding_left>=0 && ft.p3.coding_left>=0)
    {
        if((ft.p5.coding_left+ft.p3.coding_left)%3==0)
        {
            return 2;
        }
    }

    if((ft.p5.isCoding==0 || ft.p5.coding_left<0 ) && ft.p3.isCoding==1 && (ft.p3.coding_left==-1 || ft.p3.coding_left==0))//this has to be ==-1
    {
        return 1;
    }
    return 0;
}

int Annotate::getJunctionType(Gene& g, Reference & ref, Bedpe & b) {

    for (int i=0; i<fjtvec.size(); i++) {
        fjtvec[i].isInframe=getInframe(fjtvec[i]);
        bedpe_t bt=b.getBedpe(fjtvec[i].fusion_id-1);
        fjtvec[i].isCanonical=getCanonical(fjtvec[i],bt);

    }
    for (int i=0; i<fjtvec.size(); i++)
    {
        int p5_pos;
        if(fjtvec[i].p5.strand==0)
            p5_pos=fjtvec[i].p5.pos2;
        else
            p5_pos=fjtvec[i].p5.pos1;

        int p3_pos;
        if(fjtvec[i].p3.strand==0)
            p3_pos=fjtvec[i].p3.pos1;
        else
            p3_pos=fjtvec[i].p3.pos2;

        if(fjtvec[i].p5.isCoding==0)
        {
            fjtvec[i].p5Type="not-coding";
        }
        else
        {
            if(fjtvec[i].p5.coding_left==-1)
                fjtvec[i].p5Type="5p-UTR";
            else if(fjtvec[i].p5.coding_left==-2)
                fjtvec[i].p5Type="3p-UTR";
            else
                fjtvec[i].p5Type="coding";
        }

        if(fjtvec[i].p3.isCoding==0)
        {
            fjtvec[i].p3Type="not-coding";
        }
        else
        {
            if(fjtvec[i].p3.coding_left==-1)
                fjtvec[i].p3Type="5p-UTR";
            else if(fjtvec[i].p3.coding_left==-2)
                fjtvec[i].p3Type="3p-UTR";
            else
                fjtvec[i].p3Type="coding";
        }

    }



    return 0;


}



int Annotate::getJunctionPeptide(Gene& g, Reference & ref)
{
    for (int i=0; i<fjtvec.size(); i++)
    {
        int p5_pos;
        if(fjtvec[i].p5.strand==0)
            p5_pos=fjtvec[i].p5.pos2;
        else
            p5_pos=fjtvec[i].p5.pos1;
        int p3_pos;
        if(fjtvec[i].p3.strand==0)
            p3_pos=fjtvec[i].p3.pos1;
        else
            p3_pos=fjtvec[i].p3.pos2;
        vector<char> lvec=fjtvec[i].seq1;
        lvec.insert(lvec.end(),fjtvec[i].seq2.begin(),fjtvec[i].seq2.end());
        vector<char> peptide;
        int full;
        int left;
        getPeptide(fjtvec[i].seq1,lvec,fjtvec[i].p5.coding_start,peptide,full,left);
        fjtvec[i].peptide=string(peptide.begin(),peptide.end());
        /*if(peptide.size()>0)
        {
            cout<<peptide[peptide.size()-1]<<" "<<'X'<<endl;
            if(peptide[peptide.size()-1]!='X')
            {   
                cout<<fjtvec[i].peptide<<endl;
                fjtvec[i].peptide=fjtvec[i].peptide+peptide[peptide.size()-1];
                cout<<fjtvec[i].peptide<<endl;
            }
        }*/
        fjtvec[i].full=full;
        fjtvec[i].left=left;
    }

    return 0;

}


bool my_sort(const fusion_junction_t & a, const fusion_junction_t & b)
{
    if(a.p5.coding_left<b.p5.coding_left)
        return true;
    else
        return false;
}

bool my_equal(const fusion_junction_t & a, const fusion_junction_t & b)
{
    if(a.p5.coding_left==b.p5.coding_left)
        return true;
    else
        return false;
}


int Annotate::annotate(Bedpe & bedpe, Gene& g, Reference & ref)
{

    vector<fusion_junction_t> tmpvec;
    int bestInframe=0;
    int bestCanonical=0;

    vector<int> inserted;
    for(int i=0;i<bedpe.size();i++)
    {
        inserted.push_back(0);
    }
 

    int preindex=1;    
    if(fjtvec.size()>0)
        preindex=fjtvec[0].fusion_id;    

    for (int i=0; i<=fjtvec.size(); i++)
    {
        if(i==fjtvec.size() || fjtvec[i].fusion_id!=preindex)
        {
            //handle one fusion record
            sort(tmpvec.begin(),tmpvec.end(),my_sort);
            vector<fusion_junction_t>::iterator it;
            it=unique(tmpvec.begin(),tmpvec.end(),my_equal);
            tmpvec.resize(distance(tmpvec.begin(),it)); 
            //for //get strings
            string pepString="";
            string fullString="";
            string leftString="";
            if(tmpvec.size()>0)
            {
                pepString=tmpvec[0].peptide;
                fullString=to_string(static_cast<long long>(tmpvec[0].full));
                leftString=to_string(static_cast<long long>(tmpvec[0].left));
            }
            for(int x=1;x<tmpvec.size();x++)
            {
                pepString+=",";
                pepString+=tmpvec[x].peptide;

                fullString+=",";
                fullString+=to_string(static_cast<long long>(tmpvec[x].full));/////////

                leftString+=",";
                leftString+=to_string(static_cast<long long>(tmpvec[x].left));////////////       
            }
            //push_back
            bedpe.insertCol(preindex-1,to_string(static_cast<long long>(bestCanonical)));

            if(bestCanonical>0)
            {
                if(bestInframe>0)
                    bedpe.insertCol(preindex-1,string("1"));
                else
                    bedpe.insertCol(preindex-1,string("0"));
            }
            else
            {
                pepString="NA";
                fullString="NA";
                leftString="NA";
                bedpe.insertCol(preindex-1,string("0"));
            }
            bedpe.insertCol(preindex-1,pepString);
            bedpe.insertCol(preindex-1,fullString);
            bedpe.insertCol(preindex-1,leftString);
            inserted[preindex-1]=1;
            //next
            if(i<fjtvec.size())
            {
                tmpvec.clear();
                bestInframe=0;
                bestCanonical=0;

                if(fjtvec[i].peptide.size()>0)
                    tmpvec.push_back(fjtvec[i]);
                if(fjtvec[i].isInframe>bestInframe)
                    bestInframe=fjtvec[i].isInframe;
                if(fjtvec[i].isCanonical>bestCanonical)
                {
                    bestCanonical=fjtvec[i].isCanonical;
                }
            }
            preindex=fjtvec[i].fusion_id;
        }
        else
        {
            if(fjtvec[i].peptide.size()>0)
                tmpvec.push_back(fjtvec[i]);
            if(fjtvec[i].isInframe>bestInframe)
                bestInframe=fjtvec[i].isInframe;
            if(fjtvec[i].isCanonical>bestCanonical)
            {
                bestCanonical=fjtvec[i].isCanonical;
            }
        }
    }

    for(int i=0;i<inserted.size();i++)
    {
        if(inserted[i]==0)
        {
            bedpe.insertCol(i,"0");
            bedpe.insertCol(i,"0");
            bedpe.insertCol(i,"NA");
            bedpe.insertCol(i,"NA");
            bedpe.insertCol(i,"NA");
        }    
    }    

    return 0;

}

int Annotate::getCanonical(fusion_junction_t & fjt, bedpe_t & b)
{
    Bedpe bp;
    bedpe_t a;
    
    a.chr1=fjt.p5.chr;
    if(fjt.p5.strand==0)
    {
        a.start1=fjt.p5.pos2-1;
        a.end1=fjt.p5.pos2;
    }
    else
    {
        a.start1=fjt.p5.pos1-1;
        a.end1=fjt.p5.pos1;
    }

    a.chr2=fjt.p3.chr;
    if(fjt.p3.strand==0)
    {
        a.start2=fjt.p3.pos1-1;
        a.end2=fjt.p3.pos1;
    }
    else
    {
        a.start2=fjt.p3.pos2-1;
        a.end2=fjt.p3.pos2;
    }
    if(fjt.p5.strand==0)
        a.strand1="+";
    else
        a.strand1="-";
    if(fjt.p3.strand==0)
        a.strand2="+";
    else
        a.strand2="-";

//cout<<a.strand1<<"\t"<<a.start1<<"\t"<<a.end1<<"\t"<<a.strand2<<"\t"<<a.start2<<"\t"<<a.end2<<endl;

    int res=bp.isBedpeSame(a,b);
    if(res==1)
    {
        b.start1=a.start1;
        b.end1=a.end1;
        b.start2=a.start2;
        b.end2=a.end2;
    }
    return res;
     
}



