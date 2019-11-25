# ThinWire_MRIGradientCoilDesign
A gradient coil design tool for Magnetic resonance imaging (MRI)

This code demonstrates the use of thin wires to approximate a current
density. In gradient coil design for MRI usually only the z-component
of the magnetic field (Bz) is considered. Hence, only wires orthogonal to
z are used in this simulation. The Biot-Savart law is used to caclulte a 
sensitivity matrix, which is then used to calculate a current 
distribution. A regularization and an additional constraint is deployed to 
derive a ralizable coil design.
This method may be used for simple geometries. However, no generalization
for arbitrary surfaces (yet).

The code is written in MATLAB. However, no special toolboxes are used to 
make it compatible with Octave.

An explanatory introduction is given in ThinWire_Demo.m. Two basic coils
are described in the scripts Cylindrical_SingleLayer.m and 
CylindricalShielded.m.

Sebastian Littin
