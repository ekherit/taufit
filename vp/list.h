typedef struct Point_
 {
  double x,y,ey,c,m;
  struct Point_ *next;
 } tpoint,*Point;

typedef struct Tab_
 {
   Point chain;
   struct Tab_ *next;
 } tTab,*Tab;
 
       
extern Point CreatePoint(double x,double y,double ey);
extern void ExchangePoints(Point p1, Point p2);
extern void LinkPoint(Point v,Point p);
extern void AddPoint(Point v,double x,double y,double ey);
extern void PrintChain(Point v);
extern void ClearCoef(Point v);
extern void ClearMark_Chain(Point v);
extern void NegCoef_Chain(Point v);
extern void MulCoef_Chain(Point v,double k);
extern double Interpolate_Chain(Point v,double x);
extern void Mark_Chain(Point v,double x);
extern void FreeChain(Point v);

extern Tab CreateTab(Point p);
extern void AddToTab(Tab t,Point p);
extern void FreeTab(Tab t);
extern void PrintTab(Tab t);
extern double Interpolate_Tab(Tab t,double x);
extern void ClearMark_Tab(Tab t);
extern void Mark_Tab(Tab t,double x);

