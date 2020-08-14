## Graphics Vector Field Project
Project for CS 6491 Computer Graphics at Georgia Tech. This project was a two-man group effort with Steven Hillerman (also currently a GT student). Base code for 2D rendering was provided by the instructor: Jarek Rossignac.

This project showcases concepts from Computer Graphics, using a significant number of applications of geometry and operations from linear algebra.

This project is coded in Java and needs software called Processing for actual rendering. If you'd like to run the code yourself, download Processing at: https://processing.org/, download the code posted above, open "base2D.pde," and run it in Processing.

Depending on your screen resolution, you might want to change the window size (it is defined in the "base2D.pde" code on line 58.)

Key Commands:

'w': Compute and show the traces that yield the "longest" traces - the ones with that pierce the 
	most remaining traingles. The new vertices are also shown (the midpoint of a trace for a 
	particular triangle).
	<br>
'r': Retriangulate the mesh (this is the final output) - this also shows the traces
<br>
'v': toggles showing the new vertices
<br>
't': toggles showing the traces
<br>
'i': "undos" retriangualtion - this means that after hitting 'r' you can hit 'i' to go back to 
	the original mesh
	<br>
'm': this creates a bidirectional trace at the location of the mouse
<br>
'a' + leftclick: this adds a point at the location of the mouse (in order to change the mesh)
<br>
'1': this toggles showing particles which come out at the mouse's location
<br>
'2': this toggles generating the particles from the mouse's location (so hit 1 to toggle showing
	and then 2 to toggle generating)
	<br>
'3': this toggles a light background trace showing the overall field layout
<br>

Note that results (such as retriangulation) are stored and not computed at every frame.
