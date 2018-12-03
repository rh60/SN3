clear;clc;
rules=dir('*.txt');
addr='quad/';
for i=1:length(rules) 
    fn=rules(i).name;
    try
        [x,y,w] = importfile(fn);  
        nfn=strcat(addr,fn);
        dlmwrite(nfn,[(1+x)/2 (1+y)/2 w/4],'delimiter',' ','precision','%.16f');   
    catch E
        fn
    end
end    

