if(!settings.multipleView) settings.batchView=false;
settings.tex="pdflatex";
defaultfilename="p4-1";
if(settings.render < 0) settings.render=4;
settings.outformat="";
settings.inlineimage=true;
settings.embed=true;
settings.toolbar=false;
viewportmargin=(2,2);


import graph3;
import three;

size3(200);
currentprojection=orthographic(3,3,5);
currentlight=light(gray(0.4),specularfactor=3,viewport=true,
(-0.5,-0.25,0.45),(0.5,-0.5,0.5),(0.5,0.5,0.75));

int nb = 20, ns = 10;
real rb = 5.0, rs = 2.0;

triple torus(pair z) {

return ((rb + rs*cos(2*pi*z.x/ns))*cos(2*pi*z.y/nb),
(rb + rs*cos(2*pi*z.x/ns))*sin(2*pi*z.y/nb),
rs*sin(2*pi*z.x/ns));

}

surface site = scale3(0.1)*unitsphere;

for(int k1=0; k1<ns; ++k1) {
for(int k2=0; k2<nb; ++k2) {
draw(surface(torus((k1,k2))--torus((k1+1,k2))--torus((k1+1,k2+1))--torus((k1,k2+1))--cycle),
lightgray);
draw(torus((k1,k2))--torus((k1+1,k2)),Arrow3);
draw(torus((k1,k2))--torus((k1,k2+1)),Arrow3);
draw(shift(torus((k1,k2)))*site,red);
}
}
size(284.52756pt,284.52756pt,keepAspect=true);
