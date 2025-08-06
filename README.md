## Implementation steps for the 2D CFAR process.

* Nested for loop that slides the CUT across range doppler map
* For every iteration sum the signal level within all the training cells. 
* To sum convert the value from logarithmic to linear using db2pow function.
* Average the summed values for all of the training cells used. 
* After averaging convert it back to logarithimic using pow2db.
* Add the offset to it to determine the threshold. 
* Compare the signal under CUT with this threshold. If the CUT level > threshold assign
it a value of 1, else equate it to 0


## Selection of Training, Guard cells and offset

* The offset sets how far above the noise estimate the threshold should be.
    * Too low (e.g., 2–3 dB) → high false alarm rate.
    * Too high (e.g., 12+ dB) → missed detections.
* Adjustment strategy:
  * Start with 6 dB.
  * Increase by 1–2 dB increments if there are too many false alarms.
  * Decrease if real targets are being missed.

* Iterate the above process:
* Try a few combinations.
* Plot the CFAR result (overlay detected points on RDM).
* Evaluate visually and quantitatively.

## Steps taken to suppress the non-thresholded cells at the edges.

* Initialize CFAR_output as a zero matrix the same size as RDM.
* Apply CFAR thresholding only for cells where full training and guard windows exist.
* Leave edge cells (those near borders) as zeros.
* After CFAR processing, explicitly set edge regions to 0 to suppress undefined areas.