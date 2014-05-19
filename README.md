Rubik's Cube
============
Rubik's cube for [Interactive Computer Graphics](http://faculty.cooper.edu/sable2/courses/spring2013/ece462/), built with code from [Angel](http://www.cs.unm.edu/~angel/BOOK/INTERACTIVE_COMPUTER_GRAPHICS/SIXTH_EDITION/CODE/), featuring a user interface using an arcball for rotation.

setup
-----
gcc is needed with a graphics card that supports OpenGL 3. In the root directory, run

```
$ make
```

run
---

```
$ rubik.exe <filename to save to> <number of rotations | filename to load from>
```

commands
--------
Left Mouse Button - Orbit around cube
Right Mouse Button - Rotate layer
Scrollwheel - Zoom in/out

Enter - randomly rotate one layer
s - save state to the specified file
Space - reset orientation
q - quit

interface
---------
Light background - unsolved
Dark background - solved