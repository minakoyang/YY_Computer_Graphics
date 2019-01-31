# YY_Computer_Graphics
Programming Assignment: Phong Shading.

## Getting Started
Open a terminal in the folder contains the cpp, make, and model object files.
Type "make" to make the executive file.
Then type "./assn3 filename" to start the program. By default, filename is "model1.obj".
The user can load either model1 or model2 for the Model Viewer program.

The default view is the orthogonal projection with Z-buffer of the model object.
Tap 'X' to toggle between Phong Smoothing Render and Flat Shaded Render
Tap 'B' to toggle between Bump Mapping with a Dimple in the center of each triangle and Normal Mapping
Tap 'Z' to toggle between Z-buffer (colored) and wireframe view anytime.
Tap 'V' to toggle between Orthogonal mode and Perspective mode anytime.
Tap 'T' to start Translation mode:
    Tap 'A': move left;
    Tap 'D': move right;
    Tap 'W': move up;
    Tap 'S': move down.
Tap 'E' to start Scale mode:
    Tap 'A' or 'D': scale up;
    Tap 'W' or 'S': scale down.
Tap 'R' to start Rotation mode:
    Tap 'W': rotate about Y-axis counterclockwise;
    Tap 'S': rotate about Y-axis clockwise;
    Tap 'D': rotate about X-axis counterclockwise;
    Tap 'A': rotate about X-axis clockwise.


Tap 1 to reload and escape from the 'T', 'E' or 'R' mode. Z-buffer colors are randomly reassigned to the model object.
Tap esc to exit the program.

Mode information, operations, and overflow messages will show in the terminal.
