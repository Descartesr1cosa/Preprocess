#include "workflow/MPCNS_Pre_CubicSphereGridGenerator.h"
#include "workflow/MPCNS_Pre_CubicGridGenerator.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {
struct Point { double x,y,z; };
struct Block { int n, layer, face; std::string name; std::vector<Point> p; };
struct Range { int lo[3], hi[3]; };
struct Boundary { Range own; bool physical; int value; Range target; };

int IV(Param&p,const std::string&k,int v){return p.HasInt(k)?p.GetInt(k):v;}
double DV(Param&p,const std::string&k,double v){return p.HasDou(k)?p.GetDou(k):v;}
std::string SV(Param&p,const std::string&k,const std::string&v){return p.HasStr(k)?p.GetStr(k):v;}
void Fatal(const std::string&s){std::cout<<"#Fatal Error: cubic-sphere grid: "<<s<<std::endl;std::exit(-1);}
double Smooth(double s){s=std::max(0.0,std::min(1.0,s));return s*s*(3.0-2.0*s);}

Point Direction(int f,double u,double v){
    Point q;
    switch(f){
    case 0:q={1,u,v};break; case 1:q={-1,-u,v};break;
    case 2:q={-u,1,v};break; case 3:q={u,-1,v};break;
    case 4:q={-v,u,1};break; default:q={v,u,-1};break;
    }
    double d=std::sqrt(q.x*q.x+q.y*q.y+q.z*q.z); return {q.x/d,q.y/d,q.z/d};
}

std::vector<double> Geometric(double a,double b,int cells,double growth){
    if(cells<1||!(b>a)||!(growth>0)) Fatal("each radial segment needs cells>=1, end>begin, growth>0");
    std::vector<double> x(cells+1,a);
    if(std::fabs(growth-1)<1e-12){for(int i=1;i<=cells;++i)x[i]=a+(b-a)*i/cells;return x;}
    double sum=0,w=1; for(int i=0;i<cells;++i){sum+=w;w*=growth;}
    w=(b-a)/sum; for(int i=1;i<=cells;++i){x[i]=x[i-1]+w;w*=growth;} x.back()=b; return x;
}

std::vector<double> FluidCoords(Param&p){
    int count=IV(p,"cs_fluid_segments",1); if(count<1)Fatal("cs_fluid_segments must be >=1");
    std::vector<double>x(1,0); double a=0;
    for(int s=1;s<=count;++s){
        std::string id=std::to_string(s); double b=DV(p,"cs_fluid_segment_"+id+"_end",double(s)/count);
        if(b>1+1e-12)Fatal("fluid segment ends use normalized coordinates (0,1]");
        std::vector<double> part=Geometric(a,b,IV(p,"cs_fluid_segment_"+id+"_cells",16),DV(p,"cs_fluid_segment_"+id+"_growth",1.08));
        x.insert(x.end(),part.begin()+1,part.end()); a=b;
    }
    if(std::fabs(a-1)>1e-10)Fatal("last fluid segment end must be 1.0"); return x;
}

size_t Id(int n,int i,int j,int k){return static_cast<size_t>((k*n+j)*n+i);}
double PanelAngle(int index,int n){
    const double quarter=std::acos(-1.0)/4.0;const int center=n/2;
    if(index<=center)return -quarter+quarter*index/center;
    return quarter*(index-center)/(n-1-center);
}
Block SolidBlock(int face,int n,const std::vector<double>&r,const std::string&name){
    Block b={n,0,face,name,std::vector<Point>(n*n*r.size())};
    for(int k=0;k<(int)r.size();++k)for(int j=0;j<n;++j)for(int i=0;i<n;++i){
        Point d=Direction(face,std::tan(PanelAngle(i,n)),std::tan(PanelAngle(j,n)));
        b.p[Id(n,i,j,k)]={r[k]*d.x,r[k]*d.y,r[k]*d.z};} return b;
}
Block FluidBlock(int face,int n,const std::vector<double>&s,double r1,double cx,double cy,double cz,double ax,double ay,double az,const std::string&name){
    Block b={n,1,face,name,std::vector<Point>(n*n*s.size())};
    for(int k=0;k<(int)s.size();++k){double w=Smooth(s[k]),sx=r1+w*(ax-r1),sy=r1+w*(ay-r1),sz=r1+w*(az-r1);
        for(int j=0;j<n;++j)for(int i=0;i<n;++i){Point d=Direction(face,std::tan(PanelAngle(i,n)),std::tan(PanelAngle(j,n)));
            b.p[Id(n,i,j,k)]={w*cx+sx*d.x,w*cy+sy*d.y,w*cz+sz*d.z};}}
    return b;
}

Range FaceRange(int n,int nk,int f){Range r={{1,1,1},{n,n,nk}};int d=f/2;r.lo[d]=r.hi[d]=(f%2)?(d==2?nk:n):1;return r;}
Point At(const Block&b,const int q[3]){return b.p[Id(b.n,q[0]-1,q[1]-1,q[2]-1)];}
bool Same(const Point&a,const Point&b,double scale){double e=2e-11*std::max(1.0,scale);return std::fabs(a.x-b.x)<e&&std::fabs(a.y-b.y)<e&&std::fabs(a.z-b.z)<e;}
std::vector<Point> Corners(const Block&b,const Range&r,double scale){
    std::vector<Point>c; for(int m=0;m<8;++m){int q[3];for(int d=0;d<3;++d)q[d]=(m&(1<<d))?r.hi[d]:r.lo[d];Point x=At(b,q);bool dup=false;for(auto&v:c)if(Same(v,x,scale))dup=true;if(!dup)c.push_back(x);}return c;
}
bool Match(const Block&a,const Range&ra,const Block&b,const Range&rb,double scale){
    auto ca=Corners(a,ra,scale),cb=Corners(b,rb,scale);if(ca.size()!=cb.size())return false;
    for(auto&x:ca){bool found=false;for(auto&y:cb)if(Same(x,y,scale))found=true;if(!found)return false;}return true;
}
Range Target(const Block&a,const Range&ra,const Block&b,const Range&rb,double scale){
    Range out=rb;Point first=At(a,ra.lo),last=At(a,ra.hi);bool one=false,two=false;
    for(int m=0;m<8;++m){int q[3];for(int d=0;d<3;++d)q[d]=(m&(1<<d))?rb.hi[d]:rb.lo[d];Point x=At(b,q);
        if(Same(first,x,scale)){for(int d=0;d<3;++d)out.lo[d]=q[d];one=true;}
        if(Same(last,x,scale)){for(int d=0;d<3;++d)out.hi[d]=q[d];two=true;}}
    if(!one||!two)Fatal("cannot determine interface orientation");return out;
}
void MarkOrientation(const Block&a,Range&ra,const Block&b,Range&rb,double scale){
    int normal=-1;for(int d=0;d<3;++d)if(ra.lo[d]==ra.hi[d])normal=d;
    int sd=(normal==0)?1:0;
    int qa[3]={ra.lo[0],ra.lo[1],ra.lo[2]};qa[sd]=ra.hi[sd];Point probe=At(a,qa);
    int first[3]={rb.lo[0],rb.lo[1],rb.lo[2]},mapped[3]={0,0,0};bool found=false;
    Range target_face=rb;for(int d=0;d<3;++d){target_face.lo[d]=std::min(rb.lo[d],rb.hi[d]);target_face.hi[d]=std::max(rb.lo[d],rb.hi[d]);}
    for(int m=0;m<8;++m){int q[3];for(int d=0;d<3;++d)q[d]=(m&(1<<d))?target_face.hi[d]:target_face.lo[d];
        if(Same(probe,At(b,q),scale)){for(int d=0;d<3;++d)mapped[d]=q[d];found=true;}}
    if(!found)Fatal("cannot encode interface orientation");int td=-1;for(int d=0;d<3;++d)if(mapped[d]!=first[d])td=d;
    if(td<0)Fatal("degenerate interface orientation");ra.lo[sd]=-ra.lo[sd];ra.hi[sd]=-ra.hi[sd];rb.lo[td]=-rb.lo[td];rb.hi[td]=-rb.hi[td];
}

void WriteGrid(const std::string&path,const std::vector<Block>&b,const std::vector<int>&nk){
    std::ofstream f(path);if(!f)Fatal("cannot open "+path);f<<std::setprecision(16)<<b.size()<<'\n';
    for(size_t z=0;z<b.size();++z)f<<b[z].n<<' '<<b[z].n<<' '<<nk[z]<<'\n';
    for(auto&z:b)for(int c=0;c<3;++c)for(auto&q:z.p)f<<(c==0?q.x:(c==1?q.y:q.z))<<'\n';
}
void WriteInp(const std::string&path,const std::vector<Block>&blocks,const std::vector<int>&nk,bool solid,double scale){
    std::ofstream f(path);if(!f)Fatal("cannot open "+path);f<<"# multi-block equiangular cubed-sphere grid\n"<<blocks.size()<<'\n';
    for(size_t ib=0;ib<blocks.size();++ib){std::vector<Boundary>bc;
        for(int face=0;face<6;++face){Range own=FaceRange(blocks[ib].n,nk[ib],face);bool inner=face==4,outer=face==5;
            if(inner&&blocks[ib].layer==0)bc.push_back({own,true,6,{}});
            else if(inner&&blocks[ib].layer==1&&!solid)bc.push_back({own,true,6,{}});
            else if(outer&&blocks[ib].layer==1){
                const int n=blocks[ib].n,k=nk[ib]-1;
                // Classify every surface cell using its physical center. Consecutive
                // cells in an i-row are merged, so asymmetric ellipsoids need no
                // artificial logical x=0 grid line.
                for(int j=0;j<n-1;++j){
                    int begin=0,kind=0;
                    for(int i=0;i<n-1;++i){
                        const double xc=0.25*(blocks[ib].p[Id(n,i,j,k)].x+blocks[ib].p[Id(n,i+1,j,k)].x+
                                              blocks[ib].p[Id(n,i,j+1,k)].x+blocks[ib].p[Id(n,i+1,j+1,k)].x);
                        const int next=(xc>=0.0)?4:6;
                        if(i==0){kind=next;begin=0;}
                        if(next!=kind){Range patch={{begin+1,j+1,nk[ib]},{i+1,j+2,nk[ib]}};bc.push_back({patch,true,kind,{}});begin=i;kind=next;}
                        if(i==n-2){Range patch={{begin+1,j+1,nk[ib]},{n,j+2,nk[ib]}};bc.push_back({patch,true,kind,{}});}
                    }
                }
            }
            else{bool found=false;for(size_t jb=0;jb<blocks.size()&&!found;++jb)if(jb!=ib)for(int tf=0;tf<6&&!found;++tf){Range tar=FaceRange(blocks[jb].n,nk[jb],tf);
                    if(Match(blocks[ib],own,blocks[jb],tar,scale)){Range marked_own=own,marked_target=Target(blocks[ib],own,blocks[jb],tar,scale);
                        MarkOrientation(blocks[ib],marked_own,blocks[jb],marked_target,scale);bc.push_back({marked_own,false,(int)jb+1,marked_target});found=true;}}
                if(!found)Fatal("unmatched interface on block "+std::to_string(ib+1));}}
        f<<blocks[ib].n<<' '<<blocks[ib].n<<' '<<nk[ib]<<' '<<blocks[ib].name<<'\n'<<bc.size()<<'\n';
        for(auto&q:bc){for(int d=0;d<3;++d)f<<q.own.lo[d]<<' '<<q.own.hi[d]<<' ';
            if(q.physical)f<<q.value;
            else{f<<0;for(int d=0;d<3;++d)f<<' '<<q.target.lo[d]<<' '<<q.target.hi[d];f<<' '<<q.value;}
            f<<'\n';}}
}
void WriteFvbnd(const std::string&path){
    std::ofstream f(path);if(!f)Fatal("cannot open "+path);
    f<<"FVBND 1 3\n"
     <<"No_Boundary_Condition\nSolid_Surface\nSymmetry\nFarfield\nInflow\nOutflow\nPole\n"
     <<"Generic_#1\nGeneric_#2\nGeneric_#3\n"
     <<"PERIOD-X\nperiod_X\nPERIOD-Y\nperiod_Y\nPERIOD-Z\nperiod_Z\nBOUNDARIES\n";
}
}

bool CubicSphereGridGenerator::IsEnabled(Param&p){return p.HasBoo("generate_cubic_sphere_grid")&&p.GetBoo("generate_cubic_sphere_grid");}
void CubicSphereGridGenerator::GenerateIfEnabled(Param&p){
    if(!IsEnabled(p))return;int n=IV(p,"cs_surface_points",25);double r0=DV(p,"cs_r0",.5),r1=DV(p,"cs_r1",1);
    bool solid=p.HasBoo("cs_include_solid")?p.GetBoo("cs_include_solid"):true;
    const bool has_range=p.HasDou("cs_outer_xmin")||p.HasDou("cs_outer_xmax")||p.HasDou("cs_outer_ymin")||p.HasDou("cs_outer_ymax")||p.HasDou("cs_outer_zmin")||p.HasDou("cs_outer_zmax");
    double ax=DV(p,"cs_outer_axis_x",5),ay=DV(p,"cs_outer_axis_y",5),az=DV(p,"cs_outer_axis_z",5),cx=0,cy=0,cz=0;
    if(has_range){
        const double xmin=DV(p,"cs_outer_xmin",-ax),xmax=DV(p,"cs_outer_xmax",ax);
        const double ymin=DV(p,"cs_outer_ymin",-ay),ymax=DV(p,"cs_outer_ymax",ay);
        const double zmin=DV(p,"cs_outer_zmin",-az),zmax=DV(p,"cs_outer_zmax",az);
        if(!(xmin<0&&xmax>0&&ymin<ymax&&zmin<zmax))Fatal("outer bounds must satisfy xmin<0<xmax, ymin<ymax and zmin<zmax");
        cx=.5*(xmin+xmax);cy=.5*(ymin+ymax);cz=.5*(zmin+zmax);
        ax=.5*(xmax-xmin);ay=.5*(ymax-ymin);az=.5*(zmax-zmin);
    }
    if(n<3||!(r1>0)||(solid&&!(r0>0&&r1>r0)))Fatal("cs_surface_points must be >=3; radii must be valid");
    if(!(cx-ax<-r1&&cx+ax>r1&&cy-ay<-r1&&cy+ay>r1&&cz-az<-r1&&cz+az>r1))Fatal("outer ellipsoid bounds must enclose the r1 sphere");
    if(CubicPeriodicGridGenerator::IsEnabled(p))Fatal("enable only one automatic grid generator");
    std::vector<double>fluid=FluidCoords(p),sr;if(solid)sr=Geometric(r0,r1,IV(p,"cs_solid_radial_cells",12),DV(p,"cs_solid_growth",1));
    std::vector<Block>blocks;std::vector<int>nk;
    if(solid)for(int face=0;face<6;++face){blocks.push_back(SolidBlock(face,n,sr,SV(p,"cs_solid_block_name","Solid")));nk.push_back(sr.size());}
    for(int face=0;face<6;++face){blocks.push_back(FluidBlock(face,n,fluid,r1,cx,cy,cz,ax,ay,az,SV(p,"cs_fluid_block_name","Fluid")));nk.push_back(fluid.size());}
    double scale=std::max(std::fabs(cx)+ax,std::max(std::fabs(cy)+ay,std::fabs(cz)+az));WriteGrid(p.GetStr("gfilename"),blocks,nk);WriteInp(p.GetStr("bfilename"),blocks,nk,solid,scale);WriteFvbnd(p.GetStr("ffilename"));
    p.AddParam("grd_readtype",0);p.AddParam("dimension",3);if(!p.HasBoo("if_split_group_info"))p.AddParam("if_split_group_info",false);
    std::cout<<"\t-->Cubed-sphere grid generated: "<<blocks.size()<<" blocks, "<<n<<"x"<<n<<" surface points per panel.\n";
}
