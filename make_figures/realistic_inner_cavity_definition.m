%define the inner cavity in the realistic runs. Brings a1, b1, a2, b2 into global scope

a1 =  [3.3, 3.2482, 5.5487, 5.7]*1e4; b1 = [ 1.6625, 1.6361, 1.6307, 1.6625]*1e6; %north half in stereographic
a2 = [3.2482, 2.6036, 5.2548, 5.5487]*1e4; b2 = [1.6361, 1.6231, 1.6177, 1.6307]*1e6; %south half in stereographic

%single inner cavity
%a1 =  1.0e+04 *[3.5928,3.5421,   3.2584,    2.8328,    2.4983,    5.0420,    6.1061,    5.5588];
%b1 =    1.0e+06 *[    1.6648,    1.6449,    1.6337,    1.6249,    1.6222,    1.6106,    1.6235,    1.6578];
%a2 = zeros(4,1);
%b2 = zeros(4,1);


%a1 = [2.8, 2.9, 5.5487, 5.7]*1e4; b1 = [1.6625, 1.6361, 1.6307, 1.6625]*1e6; %north half
%a2 = [2.9, 2.3036, 5.2548, 5.5487]*1e4; b2 = [1.6361, 1.6231, 1.615, 1.6307]*1e6; %south half
