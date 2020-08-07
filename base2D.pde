// Template for 2D projects
// Author: Jarek ROSSIGNAC
// CS6497: Computational Aesthetics, Fall 2019, Project 3
// Students: Steven Hillerman and Michael C. Hayes
import processing.pdf.*;    // to save screen shots as PDFs, does not always work: accuracy problems, stops drawing or messes up some curves !!!
import java.awt.Toolkit;
import java.awt.datatransfer.*;

//**************************** global variables ****************************
pts P = new pts(); // class containing array of points, used to manipulate arrows
pts particles = new pts();
float t=0.5;
boolean animate=false,
        cubic=true,
        fill=false,
        genMouseParticles=false,
        improvedTraces=true,
        isRetriangulated=false,
        mouseParticles=true,
        quintic=true,
        retriangulate=true,
        showArrow=true,
        timing=false,
        showKeyArrow=true,
        spiralAverage=true,
        showFine=false,
        showCoPVectors=false,
        showMouseTrace=false,
        showMesh=true,
        showFirstField=false,
        showTrace=true,
        showTriangles=true,
        showBackgroundTraces=false,
        tracingMouse = false,
        showNewVerts = true; // toggles to display vector interpoations
int ms=0, me=0; // milli seconds start and end for timing
int npts=20000; // number of points
ARROWRING Aring = new ARROWRING();
ARROWRING RefinedAring = new ARROWRING();
ARROWRING TempAring = new ARROWRING();
int refineCounter = 6;
int f=0, df=int(pow(2,refineCounter));
float ft=0;
PFont bigFont; // for showing large labels at corner

int exitThrough=0;
MESH M = new MESH();
int cc=0; // current corner (saved, since we rebuild M at each frame)

int spawnCounter = -1;
int spawnRate = 15; // higher is slower
int maxParticles = 50;

//**************************** initialization ****************************
void setup()               // executed once at the begining 
  {
  //size(800, 800, P2D);            // window size
  size(1200, 1200, P2D);            // window size
  //frameRate(30);             // render 30 frames per second
  smooth();                  // turn on antialiasing
  P.declare(); // declares all points in P. MUST BE DONE BEFORE ADDING POINTS 
  // P.resetOnCircle(4); // sets P to have 4 points and places them in a circle on the canvas
  P.loadPts("data/pts");  // loads points form file saved with this program
  Aring.declare();
  RefinedAring.declare();
  TempAring.declare();
  //myFace = loadImage("data/pic.jpg");  // load image from file pic.jpg in folder data *** replace that file with your pic of your own face
  stevenPic = loadImage("data/steven.jpg");
  chadPic = loadImage("data/chad.png");
  textureMode(NORMAL);
  bigFont = createFont("Arial", 20); 
  textFont(bigFont); 
  textAlign(CENTER, CENTER);
  cc=M.c;
  M.reset(); 
  M.loadFromPTS(P); // loads vertices and field vectors from the sequence P of points
  M.triangulate();
  M.computeO();
  M.classifyVertices();
  M.saveData();
  
  particles.declare();
  } // end of setup

//**************************** display current frame ****************************
void draw()      // executed at each frame
  {
  if(recordingPDF) startRecordingPDF(); // starts recording graphics to make a PDF
  background(white); // clear screen and paints white background

  // ==================== MAKE ARROWS ====================   
  Aring.empty();
  for(int i=0; i<P.nv; i+=2) {Aring.addArrow(P.G[i],V(P.G[i],P.G[i+1]));}
   
  // ==================== ANIMATION ====================   
  int tm=120;
  if(animate) f=(f+1)%(tm);
  else f=(floor(ft)+tm)%(tm);
  float tt = float(f)/tm;
  t=(1.-cos(tt*TWO_PI))/2;
  
  spawnCounter = (spawnCounter + 1) % spawnRate;
 
  // ==================== TRIANGLE MESH ====================   
  
  if (showBackgroundTraces) M.showBackgroundTraces();
  
  if (showMesh)
  {
    pen(blue,2); M.drawArrows();
    noFill(); pen(black,2); M.showTriangles();
    stroke(red); M.showBorderEdges();
    M.c=cc;
    M.showCorners(3);
    noFill(); pen(black,2); M.showCurrentCorner(7);
  }
  
  if (showTrace) {
    M.showTraces();
  }
  if (showNewVerts) M.showNewVerts();
  
  // Particles from mouse
  if (mouseParticles) {
    float stepSize = 3;
    if (genMouseParticles) 
      if (spawnCounter == 0) particles.addPt(Mouse()); // Add new points
    for (int i = 0; i < particles.nv; i++) {
      int t = M.getRawContainingTriangle(particles.G[i]);
      if (t == -1) particles.removePtAt(i);
      else {
        vec P_ = M.getVectorAt(particles.G[i], t);
        pt q = P(particles.G[i], P_);
        vec q_ = M.getVectorAt(q, t);
        pt r = P(q, M(q_));
        vec S = V(particles.G[i], r);
        pt p_next = P(q, 0.5, S);
        vec dir = V(particles.G[i], p_next).normalize();
        particles.G[i] = P(particles.G[i], stepSize, dir);
      }
      if (particles.nv >= maxParticles) particles.removeFirstPt();
      fill(red);
      noStroke();
      show(particles.G[i], 6);
    }
  } //<>//
   //<>//
   // ==================== DRAW ARROWS BETWEEN CONSECUTIVE POINTS OF P ====================   
  fill(black); stroke(black);
  if(showKeyArrow) P.drawArrows(); // draws all control arrows  

  // ==================== SHOW POINTER AT MOUSE ====================   
  pt End = P(Mouse(),1,V(-2,3)), Start = P(End,20,V(-2,3)); // show semi-opaque grey arrow pointing to mouse location (useful for demos and videos)
  strokeWeight(5);  fill(grey,70); stroke(grey,70); arrow(Start,End); noFill(); 
  
  if (showMouseTrace) {
    pt mouse = P(mouseX, mouseY);
    int t = M.getContainingTriangle(mouse, 0);
    stroke(brown);
    M.traceFromPoint(mouse, t, true, false, false, false, 0);
    M.traceFromPoint(mouse, t, true, false, false, false, 0, true);
  }

  if(recordingPDF) endRecordingPDF();  // end saving a .pdf file with the image of the canvas

  fill(black); displayHeader(); // displays header
  if(scribeText && !filming) displayFooter(); // shows title, menu, and my face & name 

  if(filming && (animate || change)) snapFrameToTIF(); // saves image on canvas as movie frame 
  if(snapTIF) snapPictureToTIF();   
  if(snapJPG) snapPictureToJPG();   
  if(scribeText) {background(255,255,200); displayMenu();}
  change=false; // to avoid capturing movie frames when nothing happens
  }  // end of draw
  
