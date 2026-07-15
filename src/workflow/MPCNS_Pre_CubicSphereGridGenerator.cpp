#include "workflow/MPCNS_Pre_CubicSphereGridGenerator.h"
#include "workflow/MPCNS_Pre_CubicGridGenerator.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
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
Point Add(const Point&a,const Point&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
Point Sub(const Point&a,const Point&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
Point Mul(const Point&a,double s){return {a.x*s,a.y*s,a.z*s};}
double Dot(const Point&a,const Point&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
double Norm(const Point&a){return std::sqrt(Dot(a,a));}
Point Unit(const Point&a){const double n=Norm(a);return {a.x/n,a.y/n,a.z/n};}
Point Hermite(const Point&p0,const Point&p1,const Point&m0,const Point&m1,double t){
    const double t2=t*t,t3=t2*t;
    return Add(Add(Mul(p0,2*t3-3*t2+1),Mul(m0,t3-2*t2+t)),Add(Mul(p1,-2*t3+3*t2),Mul(m1,t3-t2)));
}

double ArcLengthParameter(const Point&p0,const Point&p1,const Point&m0,const Point&m1,double fraction){
    if(fraction<=0)return 0;if(fraction>=1)return 1;
    const int samples=256;double cumulative[samples+1];cumulative[0]=0;Point previous=p0;
    for(int q=1;q<=samples;++q){const double t=double(q)/samples;Point current=Hermite(p0,p1,m0,m1,t);cumulative[q]=cumulative[q-1]+Norm(Sub(current,previous));previous=current;}
    const double target=fraction*cumulative[samples];int hi=1;while(hi<samples&&cumulative[hi]<target)++hi;
    const double portion=(target-cumulative[hi-1])/(cumulative[hi]-cumulative[hi-1]);return (hi-1+portion)/samples;
}
double HermiteLength(const Point&p0,const Point&p1,const Point&m0,const Point&m1){
    const int samples=256;double length=0;Point previous=p0;
    for(int q=1;q<=samples;++q){Point current=Hermite(p0,p1,m0,m1,double(q)/samples);length+=Norm(Sub(current,previous));previous=current;}return length;
}
Point RayEllipsoidPoint(const Point&d,double cx,double cy,double cz,double ax,double ay,double az){
    const double aa=d.x*d.x/(ax*ax)+d.y*d.y/(ay*ay)+d.z*d.z/(az*az);
    const double bb=-2*(cx*d.x/(ax*ax)+cy*d.y/(ay*ay)+cz*d.z/(az*az));
    const double cc=cx*cx/(ax*ax)+cy*cy/(ay*ay)+cz*cz/(az*az)-1;
    const double disc=bb*bb-4*aa*cc;if(disc<=0)Fatal("a cubed-sphere ray does not intersect the outer ellipsoid");
    const double lambda=(-bb+std::sqrt(disc))/(2*aa);return Mul(d,lambda);
}
void CurveDefinition(const Point&d,double r1,double cx,double cy,double cz,double ax,double ay,double az,double tangent_scale,
                     Point&p0,Point&p1,Point&m0,Point&m1){
    p0=Mul(d,r1);p1=RayEllipsoidPoint(d,cx,cy,cz,ax,ay,az);const Point chord=Sub(p1,p0);
    const Point outer_normal=Unit({(p1.x-cx)/(ax*ax),(p1.y-cy)/(ay*ay),(p1.z-cz)/(az*az)});const double length=Norm(chord);
    const double l0=tangent_scale*std::max(0.15*length,std::min(1.5*length,Dot(chord,d)));
    const double l1=tangent_scale*std::max(0.15*length,std::min(1.5*length,Dot(chord,outer_normal)));
    m0=Mul(d,l0);m1=Mul(outer_normal,l1);
}
std::vector<double> AdvancingFractions(double length,int cells,double first_spacing,double max_ratio){
    if(cells<1||first_spacing<=0||max_ratio<1)Fatal("invalid advancing-layer controls");
    auto total=[&](double q){return std::fabs(q-1)<1e-12?first_spacing*cells:first_spacing*(std::pow(q,cells)-1)/(q-1);};
    const double effective_ratio=1.0+0.94*(max_ratio-1.0);
    double lo=1,hi=effective_ratio;
    if(total(hi)+1e-12<length)Fatal("cs_fluid_radial_cells is too small for cs_fluid_max_spacing_ratio");
    if(total(lo)>length){lo=0.01;hi=1.0;}for(int it=0;it<80;++it){const double q=.5*(lo+hi);if(total(q)<length)lo=q;else hi=q;}
    const double q=.5*(lo+hi);std::vector<double>x(cells+1,0);double width=first_spacing,sum=0;
    for(int i=0;i<cells;++i){sum+=width;x[i+1]=sum/length;width*=q;}x.back()=1;return x;
}

std::vector<double> Geometric(double a,double b,int cells,double growth){
    if(cells<1||!(b>a)||!(growth>0)) Fatal("each radial segment needs cells>=1, end>begin, growth>0");
    std::vector<double> x(cells+1,a);
    if(std::fabs(growth-1)<1e-12){for(int i=1;i<=cells;++i)x[i]=a+(b-a)*i/cells;return x;}
    double sum=0,w=1; for(int i=0;i<cells;++i){sum+=w;w*=growth;}
    w=(b-a)/sum; for(int i=1;i<=cells;++i){x[i]=x[i-1]+w;w*=growth;} x.back()=b; return x;
}

std::vector<double> FluidCoords(Param&p){
    const bool automatic=p.HasBoo("cs_fluid_auto_spacing")&&p.GetBoo("cs_fluid_auto_spacing");
    if(automatic){
        const int cells=IV(p,"cs_fluid_radial_cells",44);const double limit=DV(p,"cs_fluid_max_spacing_ratio",1.10);
        if(cells<2||limit<1.0)Fatal("automatic spacing needs cs_fluid_radial_cells>=2 and cs_fluid_max_spacing_ratio>=1");
        std::vector<double> shape(cells),width(cells),x(cells+1,0.0);const double pi=std::acos(-1.0);double max_delta=0;
        for(int i=0;i<cells;++i){const double s=std::sin(pi*(i+0.5)/cells);shape[i]=s*s;if(i)max_delta=std::max(max_delta,std::fabs(shape[i]-shape[i-1]));}
        // Leave a small allowance because the final metric uses chord lengths
        // on curved Hermite lines rather than exact parameter-space arc lengths.
        const double effective_limit=1.0+0.95*(limit-1.0);
        const double amplitude=max_delta>0?std::log(effective_limit)/max_delta:0;double sum=0;
        for(int i=0;i<cells;++i){width[i]=std::exp(amplitude*shape[i]);sum+=width[i];}
        for(int i=0;i<cells;++i)x[i+1]=x[i]+width[i]/sum;x.back()=1.0;return x;
    }
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
    // A single symmetric equiangular sequence is essential at panel seams:
    // angle(i) == -angle(n-1-i).  A piecewise sequence that forces zero into
    // an even point count breaks reversed-index seam mappings in the interior.
    const double quarter=std::acos(-1.0)/4.0;
    return -quarter+2.0*quarter*index/(n-1);
}
Block SolidBlock(int face,int n,const std::vector<double>&r,const std::string&name){
    Block b={n,0,face,name,std::vector<Point>(n*n*r.size())};
    for(int k=0;k<(int)r.size();++k)for(int j=0;j<n;++j)for(int i=0;i<n;++i){
        Point d=Direction(face,std::tan(PanelAngle(i,n)),std::tan(PanelAngle(j,n)));
        b.p[Id(n,i,j,k)]={r[k]*d.x,r[k]*d.y,r[k]*d.z};} return b;
}
Block FluidBlock(int face,int n,const std::vector<double>&s,double r1,double cx,double cy,double cz,double ax,double ay,double az,double tangent_scale,
                 bool advancing,double first_spacing,double max_ratio,const std::string&name){
    Block b={n,1,face,name,std::vector<Point>(n*n*s.size())};
    for(int j=0;j<n;++j)for(int i=0;i<n;++i){
        const Point d=Direction(face,std::tan(PanelAngle(i,n)),std::tan(PanelAngle(j,n)));Point p0,p1,m0,m1;
        CurveDefinition(d,r1,cx,cy,cz,ax,ay,az,tangent_scale,p0,p1,m0,m1);
        const std::vector<double> local=advancing?AdvancingFractions(HermiteLength(p0,p1,m0,m1),s.size()-1,first_spacing,max_ratio):s;
        for(int k=0;k<(int)local.size();++k){const double t=ArcLengthParameter(p0,p1,m0,m1,local[k]);b.p[Id(n,i,j,k)]=Hermite(p0,p1,m0,m1,t);}
    }
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

void WriteGrid(const std::string&path,const std::vector<Block>&b,const std::vector<int>&nk,int write_type){
    if(write_type<0||write_type>2)Fatal("grd_readtype must be 0 (ASCII), 1 (unformatted), or 2 (binary)");
    if(write_type==0){
        std::ofstream f(path,std::ios::out|std::ios::trunc);if(!f)Fatal("cannot open "+path);f<<std::setprecision(16)<<b.size()<<'\n';
        for(size_t z=0;z<b.size();++z)f<<b[z].n<<' '<<b[z].n<<' '<<nk[z]<<'\n';
        for(auto&z:b)for(int c=0;c<3;++c)for(auto&q:z.p)f<<(c==0?q.x:(c==1?q.y:q.z))<<'\n';
        f.close();if(!f)Fatal("failed while writing "+path);return;
    }
    std::ofstream f(path,std::ios::out|std::ios::binary|std::ios::trunc);if(!f)Fatal("cannot open "+path);
    const int count=static_cast<int>(b.size());
    if(write_type==1){const int bytes=sizeof(int);f.write(reinterpret_cast<const char*>(&bytes),sizeof(bytes));}
    f.write(reinterpret_cast<const char*>(&count),sizeof(count));
    if(write_type==1){const int bytes=sizeof(int);f.write(reinterpret_cast<const char*>(&bytes),sizeof(bytes));const int dims_bytes=3*count*sizeof(int);f.write(reinterpret_cast<const char*>(&dims_bytes),sizeof(dims_bytes));}
    for(size_t z=0;z<b.size();++z){const int dims[3]={b[z].n,b[z].n,nk[z]};f.write(reinterpret_cast<const char*>(dims),sizeof(dims));}
    if(write_type==1){const int dims_bytes=3*count*sizeof(int);f.write(reinterpret_cast<const char*>(&dims_bytes),sizeof(dims_bytes));}
    for(size_t z=0;z<b.size();++z){
        const int data_bytes=static_cast<int>(3*b[z].p.size()*sizeof(double));
        if(write_type==1)f.write(reinterpret_cast<const char*>(&data_bytes),sizeof(data_bytes));
        for(int c=0;c<3;++c)for(auto&q:b[z].p){const double value=(c==0?q.x:(c==1?q.y:q.z));f.write(reinterpret_cast<const char*>(&value),sizeof(value));}
        if(write_type==1)f.write(reinterpret_cast<const char*>(&data_bytes),sizeof(data_bytes));
    }
    f.close();if(!f)Fatal("failed while writing "+path);
}
void WriteInp(const std::string&path,const std::vector<Block>&blocks,const std::vector<int>&nk,bool solid,double scale){
    const std::string temporary=path+".tmp";
    std::ofstream f(temporary,std::ios::out|std::ios::trunc);if(!f)Fatal("cannot open "+temporary);f<<"# multi-block equiangular cubed-sphere grid\n"<<blocks.size()<<'\n';
    for(size_t ib=0;ib<blocks.size();++ib){std::vector<Boundary>bc;
        for(int face=0;face<6;++face){Range own=FaceRange(blocks[ib].n,nk[ib],face);bool inner=face==4,outer=face==5;
            if(inner&&blocks[ib].layer==0)bc.push_back({own,true,6,{}});
            else if(inner&&blocks[ib].layer==1&&!solid)bc.push_back({own,true,2,{}});
            else if(outer&&blocks[ib].layer==1){
                const int n=blocks[ib].n,k=nk[ib]-1;
                std::vector<double> cell_x((n-1)*(n-1));int positive=0,negative=0;
                for(int j=0;j<n-1;++j)for(int i=0;i<n-1;++i){
                    const double x=.25*(blocks[ib].p[Id(n,i,j,k)].x+blocks[ib].p[Id(n,i+1,j,k)].x+
                                         blocks[ib].p[Id(n,i,j+1,k)].x+blocks[ib].p[Id(n,i+1,j+1,k)].x);
                    cell_x[j*(n-1)+i]=x;if(x>=0)++positive;else ++negative;
                }
                if(negative==0)bc.push_back({own,true,4,{}});
                else if(positive==0)bc.push_back({own,true,6,{}});
                else{
                    // Find the best two-rectangle classification of surface
                    // cells. No node is required to lie on physical x=0.
                    double best=std::numeric_limits<double>::max();int split_dim=0,cut=1;bool low_farfield=true;
                    for(int dim=0;dim<2;++dim)for(int c=1;c<n-1;++c)for(int orientation=0;orientation<2;++orientation){
                        double cost=0;for(int j=0;j<n-1;++j)for(int i=0;i<n-1;++i){
                            const bool low=(dim==0?i:j)<c;const bool predicted_far=low?(orientation==0):(orientation!=0);
                            const double x=cell_x[j*(n-1)+i];if(predicted_far!=(x>=0))cost+=std::fabs(x);
                        }
                        if(cost<best){best=cost;split_dim=dim;cut=c;low_farfield=orientation==0;}
                    }
                    Range low=own,high=own;low.hi[split_dim]=cut+1;high.lo[split_dim]=cut+1;
                    bc.push_back({low,true,low_farfield?4:6,{}});
                    bc.push_back({high,true,low_farfield?6:4,{}});
                }
            }
            else{bool found=false;for(size_t jb=0;jb<blocks.size()&&!found;++jb)if(jb!=ib)for(int tf=0;tf<6&&!found;++tf){Range tar=FaceRange(blocks[jb].n,nk[jb],tf);
                    if(Match(blocks[ib],own,blocks[jb],tar,scale)){Range marked_own=own,marked_target=Target(blocks[ib],own,blocks[jb],tar,scale);
                        MarkOrientation(blocks[ib],marked_own,blocks[jb],marked_target,scale);bc.push_back({marked_own,false,(int)jb+1,marked_target});found=true;}}
                if(!found)Fatal("unmatched interface on block "+std::to_string(ib+1));}}
        f<<blocks[ib].n<<' '<<blocks[ib].n<<' '<<nk[ib]<<' '<<blocks[ib].name<<'\n'<<bc.size()<<'\n';
        for(auto&q:bc){for(int d=0;d<3;++d)f<<q.own.lo[d]<<' '<<q.own.hi[d]<<' ';
            if(q.physical)f<<q.value;
            else{f<<-1<<'\n';for(int d=0;d<3;++d)f<<q.target.lo[d]<<' '<<q.target.hi[d]<<' ';f<<q.value;}
            f<<'\n';}}
    f.close();if(!f)Fatal("failed while writing "+temporary);
    if(std::rename(temporary.c_str(),path.c_str())!=0)Fatal("cannot replace "+path+" with completed inp file");
}
void WriteFvbnd(const std::string&path){
    std::ofstream f(path,std::ios::out|std::ios::trunc);if(!f)Fatal("cannot open "+path);
    f<<"FVBND 1 3\n"
     <<"No_Boundary_Condition\nSolid_Surface\nSymmetry\nFarfield\nInflow\nOutflow\nPole\n"
     <<"Generic_#1\nGeneric_#2\nGeneric_#3\n"
     <<"PERIOD-X\nperiod_X\nPERIOD-Y\nperiod_Y\nPERIOD-Z\nperiod_Z\nBOUNDARIES\n";
    f.close();if(!f)Fatal("failed while writing "+path);
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
    std::vector<double>sr;
    if(solid){
        const double solid_growth=DV(p,"cs_solid_growth",1.0);
        if(!(solid_growth>0.0))Fatal("cs_solid_growth must be > 0");
        // User-facing growth is the inward coarsening ratio: values > 1
        // make cells progressively thinner from r0 toward the interface r1.
        sr=Geometric(r0,r1,IV(p,"cs_solid_radial_cells",12),1.0/solid_growth);
    }
    std::vector<Block>blocks;std::vector<int>nk;
    if(solid)for(int face=0;face<6;++face){blocks.push_back(SolidBlock(face,n,sr,SV(p,"cs_solid_block_name","Solid")));nk.push_back(sr.size());}
    const double tangent_scale=DV(p,"cs_orthogonal_tangent_scale",1.0);
    if(!(tangent_scale>0.0&&tangent_scale<=2.0))Fatal("cs_orthogonal_tangent_scale must be in (0,2]");
    const bool advancing=p.HasBoo("cs_fluid_advance_layers")&&p.GetBoo("cs_fluid_advance_layers");
    const double interface_ratio=DV(p,"cs_interface_spacing_ratio",1.0),max_ratio=DV(p,"cs_fluid_max_spacing_ratio",1.10);
    const double first_spacing=solid?(sr.back()-sr[sr.size()-2])*interface_ratio:DV(p,"cs_fluid_first_spacing",0.02*r1);
    std::vector<double>fluid;
    if(advancing){
        double longest=0;for(int face=0;face<6;++face)for(int j=0;j<n;++j)for(int i=0;i<n;++i){
            const Point d=Direction(face,std::tan(PanelAngle(i,n)),std::tan(PanelAngle(j,n)));Point p0,p1,m0,m1;
            CurveDefinition(d,r1,cx,cy,cz,ax,ay,az,tangent_scale,p0,p1,m0,m1);longest=std::max(longest,HermiteLength(p0,p1,m0,m1));}
        int cells=IV(p,"cs_fluid_radial_cells",0);
        if(cells==0){const double effective_ratio=1.0+0.94*(max_ratio-1.0);cells=1;double covered=first_spacing;
            while(covered<longest){++cells;covered=std::fabs(effective_ratio-1)<1e-12?first_spacing*cells:first_spacing*(std::pow(effective_ratio,cells)-1)/(effective_ratio-1);}}
        fluid.resize(cells+1);for(int k=0;k<=cells;++k)fluid[k]=double(k)/cells;
        std::cout<<"\t   advancing layers: cells="<<cells<<", first spacing="<<first_spacing<<", max growth="<<max_ratio<<".\n";
    }else fluid=FluidCoords(p);
    for(int face=0;face<6;++face){blocks.push_back(FluidBlock(face,n,fluid,r1,cx,cy,cz,ax,ay,az,tangent_scale,advancing,first_spacing,max_ratio,SV(p,"cs_fluid_block_name","Fluid")));nk.push_back(fluid.size());}
    const int grid_type=IV(p,"grd_readtype",0);
    double scale=std::max(std::fabs(cx)+ax,std::max(std::fabs(cy)+ay,std::fabs(cz)+az));WriteGrid(p.GetStr("gfilename"),blocks,nk,grid_type);WriteInp(p.GetStr("bfilename"),blocks,nk,solid,scale);WriteFvbnd(p.GetStr("ffilename"));
    p.AddParam("dimension",3);if(!p.HasBoo("if_split_group_info"))p.AddParam("if_split_group_info",false);
    std::cout<<"\t-->Cubed-sphere grid generated: "<<blocks.size()<<" blocks, "<<n<<"x"<<n<<" surface points per panel, grd_readtype="<<grid_type<<".\n";
}
