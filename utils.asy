import graph;
import palette;

//http://www.sandia.gov/~kmorel/documents/ColorMaps/CoolWarmUChar33.csv
pen[] paraview_cooltowarm=Gradient(
		       rgb(59,76,192),
		       rgb(68,90,204),
		       rgb(77,104,215),
		       rgb(87,117,225),
		       rgb(98,130,234),
		       rgb(108,142,241),
		       rgb(119,154,247),
		       rgb(130,165,251),
		       rgb(141,176,254),
		       rgb(152,185,255),
		       rgb(163,194,255),
		       rgb(174,201,253),
		       rgb(184,208,249),
		       rgb(194,213,244),
		       rgb(204,217,238),
		       rgb(213,219,230),
		       rgb(221,221,221),
		       rgb(229,216,209),
		       rgb(236,211,197),
		       rgb(241,204,185),
		       rgb(245,196,173),
		       rgb(247,187,160),
		       rgb(247,177,148),
		       rgb(247,166,135),
		       rgb(244,154,123),
		       rgb(241,141,111),
		       rgb(236,127,99),
		       rgb(229,112,88),
		       rgb(222,96,77),
		       rgb(213,80,66),
		       rgb(203,62,56),
		       rgb(192,40,47),
		       rgb(180,4,38)	       );


// Find the comma-separated strings to use in the legend
string[] set_legends(string runlegs)
{
  string[] legends;
  bool myleg=((runlegs== "") ? false: true);
  bool flag=true;
  int n=-1;
  int lastpos=0;
  string legends[];
  if(myleg) {
    string runleg;
    while(flag) {
      ++n;
      int pos=find(runlegs,",",lastpos);
      if(lastpos == -1) {runleg=""; flag=false;}
    
      runleg=substr(runlegs,lastpos,pos-lastpos);

      lastpos=pos > 0 ? pos+1 : -1;
      if(flag) legends.push(runleg);
    }
  }
  return legends;
}

// optionally draw different data from files as well:
void draw_another(bool myleg,string[] legends, int n)
{
  bool anotherlegend=true;
  bool anotherone=(getstring("compare with other file? y/n") == "y");
  while(anotherone) {
    string filename;
    filename=getstring("filename:");
    file fin=input(filename).line();
    real[][] a=fin.dimension(0,0);
    a=transpose(a);
    int n0=getint("column x");
    int n1=getint("column y");
    string legend=myleg ? legends[n] : texify(filename);
    if(anotherlegend) {
      draw(graph(a[n0],a[n1]),Pen(n),legend);
      anotherlegend=false;
    } else {
      draw(graph(a[n0],a[n1]),Pen(n));
    }
    anotherone=(getstring("another one? y/n") == "y");
  }
}


// Load file (note format change from FORTRAN to C)
real[][][] readfile(int nx, int ny, int nz, string name) {
  file fin=input(name,mode="binary").singlereal();
  //file fin=input(name,mode="binary");
  //file fin=input(name);
  real[][][] f=new real[nx][ny][nz];

  // There is some really weird ordering in .h5 files.
  for(int k=0; k < nz; ++k) {
    for(int j=0; j < ny; ++j) {
      for(int i=0; i < nx; ++i) {
	f[i][j][k]=fin.read(1);
      }
    }
  }
  return f;
}

// return a 2D cross-section of a 3D file in index dir.
real[][] cut2(real[][][] f,
    int nx, int ny, int nz, int c, int dir) {
  real[][] f2;
  if(dir == 1) {
    f2=new real[nx][ny];
    for(int i=0; i < nx; ++i) {
      for(int j=0; j < ny; ++j) {
	f2[i][j]=f[i][j][c];
      }
    }
  }
  if(dir == 2) {
    f2=new real[ny][nz];
    for(int j=0; j < ny; ++j) {
      for(int k=0; k < nz; ++k) {
	f2[j][k]=f[c][j][k];
      }
    }
  }
  if(dir == 3) {
    f2=new real[nx][nz];
    for(int i=0; i < nx; ++i) {
      for(int k=0; k < nz; ++k) {
	f2[i][k]=f[i][c][k];
      }
    }
  }
  return f2;
}

// return a 1D cross-section of a 3D file in index dir.
real[] cut1(real[][][] f, int nx, int ny, int nz, int c1, int c2, int dir) {
  real [] f1;
  if(dir == 1) {
    f1=new real[nx];
    for(int i=0; i < nx; ++i)
      f1[i]=f[i][c1][c2];
  }
  if(dir == 2) {
    f1=new real[ny];
  for(int j=0; j < ny; ++j)
    f1[j]=f[c1][j][c2];
  }
  if(dir == 3) {
    f1=new real[nz];
    for(int k=0; k < ny; ++k)
      f1[k]=f[c1][c2][k];
  }
  return f1;
}


// Given a function y(x), find the fit to a+by for x in (start,stop).
// The data is returned in a pair (a,b).
pair linear(real[] x, real[] y, real start, real stop)
{
  pair AB=(0,0);

  real meanx=0.0, meany=0.0;
  int N=0;
  for(int i=0; i < x.length; ++i) {
    if(x[i] > start && x[i] < stop) {
      meanx += x[i];
      meany += y[i];
      ++N;
    }
  }
  if(N == 0) return AB;
  meanx /= N;
  meany /= N;

  real a=0, b=0;
  for(int i=0; i < x.length; ++i) {
    if(x[i] > start && x[i] < stop) {
      a += (x[i]-meanx)*(y[i]-meany);
      b += (x[i]-meanx)*(x[i]-meanx);
    }
  }
  
  return (meany-(a/b)*meanx,a/b);
 }

void maskit(real[][][] f,int nx, int ny, int nz)
{
  string maskname=getstring("mask filename");
  real[][][] mask=readfile(nx,ny,nz,maskname);
  for(int i=0; i < nx; ++i) {
    for(int j=0; j < ny; ++j) {
      for(int k=0; k < nz; ++k) {
	if(mask[i][j][k] != 0.0)
	  f[i][j][k] = 0.0;
      }
    }
  }
}

void clipellipse(real l1, real l2, picture pic)
{
  pair O=(l1/2,l2/2);
  real ay=getreal("ay");
  if(ay != 0.0) {
    path wall=ellipse(O,1,1/sqrt(ay));
    draw(pic,wall);
    clip(pic,wall);
  }
}

void cliprectangle(real l1, real l2, picture pic)
{
  real w=getreal("width");
  real h=getreal("height");

  pair O=(l1/2,l2/2);
  pair pw=(0.5*w,0), ph=(0,0.5*h);
  path wall=O-pw-ph--O-pw+ph--O+pw+ph--O+pw-ph--cycle;
  draw(pic,wall);
  clip(pic,wall);
}

pair imagedims(string cutdir)
{
  // Get dimensions of image:
  string sl1, sl2;
  if(cutdir == "x") { sl1="yl"; sl2="zl";}
  if(cutdir == "y") { sl1="xl"; sl2="zl";}
  if(cutdir == "z") { sl1="xl"; sl2="yl";} // CHECKME
  real l1=getreal(sl1);
  real l2=getreal(sl2);
  return(l1,l2);
}