# ThinWire_MRIGradientCoilDesign
A gradient coil design tool for Magnetic resonance imaging (MRI)

This code demonstrates the use of thin wires to approximate a current
density. In gradient coil design for MRI usually only the the z-component
of the magnetic field (Bz) is considered. Hence, only wires orthogonal to
z are used in this simulation. A sensitivity matrix is used to calculate
acurrent distribution. A regularization and an additional constraint is
deployed to derive a ralizable coil design.
This ,ethod may be used for simple geometries. However, no generalization
for arbitrary surfaces.

Sebastian Littin
