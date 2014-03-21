Summary: In this project, you have to define a path for the planes appearing on the screen and make them land without colliding them.  Basically we select the plane by clicking the mouse on it and then the mouse defines the path for the plane. We have to make it land at the runway and have to make sure that planes don’t crash. If the path turns green, then the plane is set to land. 

Usage: Please run the file using runme.m

Options:

Features: The game starts with a defined number of planes.
Planes appear with a certain randomness on the screen. 
If the planes get close to each other a sound plays warning the player.
A crashing sound plays if the plane collides.

Theory:
 The plane follows the path exactly with a constant velocity. It does not depend on the speed with which we draw the path. The plane speed is independent of any other factor. This was accomplished by calculating the direction vector of next point on the path, obtaining the distance between two points and calculating the time the plane has to travel in that direction using the formula Time=Distance/Velocity. Though it is still not time dependent because the Pause may vary, causing the plane to detract. Rather, the code calculates the number of iterations it has to follow before it moves on to the next point. Also, in order to make sure that I don’t get NaN direction vectors, I made sure that the code did not plot a new point at the same place. This also prevents the game from having too much coordinates.
