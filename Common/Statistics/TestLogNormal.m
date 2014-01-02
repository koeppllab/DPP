
clear all; 
close all;

lMu = [log(3); log(30)];
lSigma = [0.5, 0.4
          0.4, 0.8];

X = DrawLogNormal(lMu + diag(lSigma), lSigma, 500000);

[h, b] = hist(X(1, :), 100);
dB = b(2)-b(1);
h = h/sum(h*dB);

for l=1:length(b)
    lPDF(l) = LogNormalPDF(lMu(1), lSigma(1), b(l));
end
lPDF = lPDF / sum(dB*lPDF);

plot(b, h, 'r'); hold on;
plot(b, lPDF, 'k');

bOld = b;


figure;



subplot(1,2,1);


[bins, centers] = hist3(X([2, 1], :)', [20, 20]);
xVec = centers{1};
yVec = centers{2};
db1 = xVec(2) - xVec(1);
db2 = yVec(2) - yVec(1);
bins = log(bins / (db1*db2*sum(sum(bins))));

lins = linspace(-30, 0, 30);

contourf(centers{2}, centers{1}, bins, lins);
colormap = 1 - colormap(gray);


for l=1:length(xVec)
    for k=1:length(yVec)
        lPDFJ(l, k) = LogNormalPDF(lMu, lSigma, [yVec(k); xVec(l);]);
    end
end
lPDFJ = lPDFJ / (db1*db2*sum(sum(lPDFJ)));
subplot(1,2,2);
contourf(centers{2}, centers{1}, log(lPDFJ), lins);



vChain = zeros(2, 100000);
vChain(:, 1) = [3, 30];

sigmaProp = 2;

for k=1:length(vChain)
   
    vProp(1) = lognrnd(log(vChain(1, k)), sigmaProp);
    vProp(2) = lognrnd(log(vChain(2, k)), sigmaProp);
    
    forwardProb = 0;
    backwardProb = 0;
    
    for i=1:length(vProp)
       forwardProb = forwardProb + log(lognpdf(vProp(i), log(vChain(i, k)), sigmaProp)); 
       backwardProb = backwardProb + log(lognpdf(vChain(i, k), log(vProp(i)), sigmaProp));
    end
    
    LNew = log(LogNormalPDF(lMu, lSigma, vProp'));
    
    if (k==1)
        a = 1;
        
    else
        
       a = min(1, exp(LNew-LOld + backwardProb - forwardProb));
      
    end
    
    if (rand<a)
       
        fprintf('accept\n');
        LOld = LNew;
        vChain(:, k+1) = vProp;
        
    else
       fprintf('reject\n');
       vChain(:, k+1) = vChain(:, k);
        
    end
    
    if (mod(k, 10000)==0)
    plot(vChain(:, 1:k)'); drawnow;
    end
end


X = DrawLogNormal(lMu, lSigma, 1000000);

[h, b] = hist(X(2, :), 100);
h = h/(sum(h)*(b(2)-b(1)));

plot(b, h); hold on;

[h, b] = hist(vChain(2, :), 100);
h = h/(sum(h)*(b(2)-b(1)));

plot(b, h, 'r'); hold on;



%X = DrawLogNormal(lMu, lSigma, 1000000);

