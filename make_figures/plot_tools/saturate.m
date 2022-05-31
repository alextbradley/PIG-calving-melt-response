function M2 = saturate(M1, levt, levb)
M2 = M1;
M2(M2 > levt) = levt;
M2(M2 < levb) = levb;
M2(1,1) = levt;
M2(2,1) = levb;

