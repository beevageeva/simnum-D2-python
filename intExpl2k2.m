$Assumptions = {Element[{k0,k1, z0, z1 , zi0, zi1, zc0, zc1,W}, Reals],  zi0 >0, zi1>0 zc0 >0, zc1>0,  k0>0, k1>0 , W>0}

Print[Integrate[ Exp[-((z0-zc0)^2-(z1-zc1)^2)/W^2]  Cos[k0 (z0 - zi0) + k1 (z1-zi1) ]  Exp[I m z0 ] Exp[I n z1],  {z0, -Infinity, Infinity},  {z1, -Infinity, Infinity} ]]

Print[FullSimplify[Integrate[Exp[-((z0-zc0)^2-(z1-zc1)^2)/W^2]  Cos[k0 (z0 - zi0) + k1 (z1-zi1) ] Exp[I m z0] Exp[I n z1],  {z0, -Infinity, Infinity}, {z1, -Infinity, Infinity} ]]]


