function ElementData=CirPlate2(n)
% 4三棱柱单元=圆板 
% 模型网格录入程序
%% element1 Prism
ele1=struct;
ele1.Type=2;
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,[n-1,n-1],[n-1,n-1],[n-1,n-1],n-2};
number.H=[n-2,n-1];
ele1.Number=number;
number=struct;
number.edge=[0,18,0,0,0,0,0,18,0];
number.face={0,[0,0],[0,0],[0,0],0};
number.H=[0,0];
ele1.MapNumber=number;
ele1.MapRight={[0;100;6.12323399573677e-15;0;100;6.12323399573677e-15;8.25793454723324;16.4594590280734;24.5485487140799;32.4699469204684;40.1695424652969;47.5947393037074;54.6948158122427;61.4212712689668;67.7281571625741;73.5723910673132;78.9140509396394;83.7166478262529;87.9473751206489;91.5773326655057;94.5817241700635;96.9400265939330;98.6361303402722;99.6584493006670;8.25793454723324;16.4594590280734;24.5485487140799;32.4699469204684;40.1695424652969;47.5947393037074;54.6948158122427;61.4212712689668;67.7281571625741;73.5723910673132;78.9140509396394;83.7166478262529;87.9473751206489;91.5773326655057;94.5817241700635;96.9400265939330;98.6361303402722;99.6584493006670],...
[0;0;100;0;0;100;99.6584493006670;98.6361303402722;96.9400265939330;94.5817241700635;91.5773326655058;87.9473751206489;83.7166478262529;78.9140509396394;73.5723910673132;67.7281571625741;61.4212712689668;54.6948158122427;47.5947393037074;40.1695424652969;32.4699469204683;24.5485487140799;16.4594590280734;8.25793454723323;99.6584493006670;98.6361303402722;96.9400265939330;94.5817241700635;91.5773326655058;87.9473751206489;83.7166478262529;78.9140509396394;73.5723910673132;67.7281571625741;61.4212712689668;54.6948158122427;47.5947393037074;40.1695424652969;32.4699469204683;24.5485487140799;16.4594590280734;8.25793454723323],...
[0;0;0;10;10;10;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10]};
%% element2 Prism
ele2=struct;
ele2.Type=2;
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,[n-1,n-1],[n-1,n-1],[n-1,n-1],n-2};
number.H=[n-2,n-1];
ele2.Number=number;
number=struct;
number.edge=[0,18,0,0,0,0,0,18,0];
number.face={0,[0,0],[0,0],[0,0],0};
number.H=[0,0];
ele2.MapNumber=number;
ele2.MapRight={[0;6.12323399573677e-15;-100;0;6.12323399573677e-15;-100;-99.6584493006670;-98.6361303402722;-96.9400265939330;-94.5817241700635;-91.5773326655058;-87.9473751206489;-83.7166478262529;-78.9140509396394;-73.5723910673132;-67.7281571625741;-61.4212712689668;-54.6948158122427;-47.5947393037073;-40.1695424652969;-32.4699469204683;-24.5485487140799;-16.4594590280734;-8.25793454723323;-99.6584493006670;-98.6361303402722;-96.9400265939330;-94.5817241700635;-91.5773326655058;-87.9473751206489;-83.7166478262529;-78.9140509396394;-73.5723910673132;-67.7281571625741;-61.4212712689668;-54.6948158122427;-47.5947393037073;-40.1695424652969;-32.4699469204683;-24.5485487140799;-16.4594590280734;-8.25793454723323],...
 [0;100;1.22464679914735e-14;0;100;1.22464679914735e-14;8.25793454723322;16.4594590280734;24.5485487140800;32.4699469204684;40.1695424652969;47.5947393037074;54.6948158122427;61.4212712689668;67.7281571625741;73.5723910673132;78.9140509396394;83.7166478262529;87.9473751206489;91.5773326655058;94.5817241700635;96.9400265939330;98.6361303402723;99.6584493006670;8.25793454723322;16.4594590280734;24.5485487140800;32.4699469204684;40.1695424652969;47.5947393037074;54.6948158122427;61.4212712689668;67.7281571625741;73.5723910673132;78.9140509396394;83.7166478262529;87.9473751206489;91.5773326655058;94.5817241700635;96.9400265939330;98.6361303402723;99.6584493006670],...
[0;0;0;10;10;10;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10]};
%% element3 Prism
ele3=struct;
ele3.Type=2;
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,[n-1,n-1],[n-1,n-1],[n-1,n-1],n-2};
number.H=[n-2,n-1];
ele3.Number=number;
number=struct;
number.edge=[0,18,0,0,0,0,0,18,0];
number.face={0,[0,0],[0,0],[0,0],0};
number.H=[0,0];
ele3.MapNumber=number;
ele3.MapRight={[0;-100;-1.83697019872103e-14;0;-100;-1.83697019872103e-14;-8.25793454723327;-16.4594590280734;-24.5485487140799;-32.4699469204684;-40.1695424652970;-47.5947393037074;-54.6948158122428;-61.4212712689668;-67.7281571625741;-73.5723910673132;-78.9140509396394;-83.7166478262529;-87.9473751206489;-91.5773326655058;-94.5817241700635;-96.9400265939330;-98.6361303402723;-99.6584493006670;-8.25793454723327;-16.4594590280734;-24.5485487140799;-32.4699469204684;-40.1695424652970;-47.5947393037074;-54.6948158122428;-61.4212712689668;-67.7281571625741;-73.5723910673132;-78.9140509396394;-83.7166478262529;-87.9473751206489;-91.5773326655058;-94.5817241700635;-96.9400265939330;-98.6361303402723;-99.6584493006670],...
[0;1.22464679914735e-14;-100;0;1.22464679914735e-14;-100;-99.6584493006670;-98.6361303402723;-96.9400265939330;-94.5817241700635;-91.5773326655057;-87.9473751206489;-83.7166478262528;-78.9140509396394;-73.5723910673132;-67.7281571625741;-61.4212712689668;-54.6948158122427;-47.5947393037074;-40.1695424652969;-32.4699469204683;-24.5485487140799;-16.4594590280734;-8.25793454723316;-99.6584493006670;-98.6361303402723;-96.9400265939330;-94.5817241700635;-91.5773326655057;-87.9473751206489;-83.7166478262528;-78.9140509396394;-73.5723910673132;-67.7281571625741;-61.4212712689668;-54.6948158122427;-47.5947393037074;-40.1695424652969;-32.4699469204683;-24.5485487140799;-16.4594590280734;-8.25793454723316],...
[0;0;0;10;10;10;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10]};
%% element4 Prism
ele4=struct;
ele4.Type=2;
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,[n-1,n-1],[n-1,n-1],[n-1,n-1],n-2};
number.H=[n-2,n-1];
ele4.Number=number;
number=struct;
number.edge=[0,18,0,0,0,0,0,18,0];
number.face={0,[0,0],[0,0],[0,0],0};
number.H=[0,0];
ele4.MapNumber=number;
ele4.MapRight={[0;-1.83697019872103e-14;100;0;-1.83697019872103e-14;100;99.6584493006670;98.6361303402723;96.9400265939330;94.5817241700635;91.5773326655057;87.9473751206489;83.7166478262529;78.9140509396393;73.5723910673131;67.7281571625741;61.4212712689668;54.6948158122427;47.5947393037074;40.1695424652970;32.4699469204683;24.5485487140799;16.4594590280733;8.25793454723324;99.6584493006670;98.6361303402723;96.9400265939330;94.5817241700635;91.5773326655057;87.9473751206489;83.7166478262529;78.9140509396393;73.5723910673131;67.7281571625741;61.4212712689668;54.6948158122427;47.5947393037074;40.1695424652970;32.4699469204683;24.5485487140799;16.4594590280733;8.25793454723324],...
[0;-100;-2.44929359829471e-14;0;-100;-2.44929359829471e-14;-8.25793454723328;-16.4594590280734;-24.5485487140799;-32.4699469204684;-40.1695424652970;-47.5947393037074;-54.6948158122427;-61.4212712689668;-67.7281571625742;-73.5723910673132;-78.9140509396394;-83.7166478262529;-87.9473751206489;-91.5773326655057;-94.5817241700635;-96.9400265939331;-98.6361303402723;-99.6584493006670;-8.25793454723328;-16.4594590280734;-24.5485487140799;-32.4699469204684;-40.1695424652970;-47.5947393037074;-54.6948158122427;-61.4212712689668;-67.7281571625742;-73.5723910673132;-78.9140509396394;-83.7166478262529;-87.9473751206489;-91.5773326655057;-94.5817241700635;-96.9400265939331;-98.6361303402723;-99.6584493006670],...
[0;0;0;10;10;10;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;10]};
           


ElementData={ele1,ele2,ele3,ele4};
% ElementData={ele1,ele2,ele3,ele4};
% ElementData={ele5,ele6,ele7,ele8};
% ElementData={ele1};

end