%-------------------------------------------------------------
%    This is the file Aggregation_Pi.m
%
%    Version Nov 2019.
%    Simone Coniglio <simone.coniglio@airbus.com>
%    Propulsion Airframe Stress Transverse,
%    31300 Toulouse, France.
%
function [Wa,dWa]=Aggregation_Pi(z,p)
% function that make the aggregation of the values in z and also compute
% sensitivities
zm=repmat(max(z),size(z,1),1);
ka=p.ka;
switch p.aggregation
    case 'asymptotic'
        Wa=sum(z);
        dWa=ones(size(z));
    case 'boolean'
        Wa=1-prod(1-z);
        dWa=repmat(prod(1-z.*(z~=1)),size(z,1),1)./(1-z.*(z~=1));
    case 'p-norm'
        zp=p.zp;
        zm=zm+zp;
        z=z+zp;
        Wa=exp(zm(1,:)).*(sum((z./exp(zm)).^ka,1)).^(1/ka)-zp;
        dWa=(z./exp(zm)).^(ka-1).*repmat((sum((z./exp(zm)).^ka,1)).^(1/ka-1),size(z,1),1);
    case 'p-mean'
        zp=p.zp;
        zm=zm+zp;
        z=z+zp;
        Wa=exp(zm(1,:)).*(mean((z./exp(zm)).^ka,1)).^(1/ka)-zp;
        dWa=1/size(z,1)^(1/ka)*(z./exp(zm)).^(ka-1).*repmat((sum((z./exp(zm)).^ka,1)).^(1/ka-1),size(z,1),1);
    case 'KS'
        Wa=zm(1,:)+1/ka*log(sum(exp(ka*(z-zm)),1));
        dWa=exp(ka*(z-zm))./repmat(sum(exp(ka*(z-zm)),1),size(z,1),1);
    case 'KSl'
        Wa=zm(1,:)+1/ka*log(mean(exp(ka*(z-zm)),1));
        dWa=exp(ka*(z-zm))./repmat(sum(exp(ka*(z-zm)),1),size(z,1),1);
    case 'IE'
        Wa=sum(z.*exp(ka*(z-zm)))./sum(exp(ka*(z-zm)),1);
        dWa=((exp(ka*(z-zm))+ka*z.*exp(ka*(z-zm))).*repmat(sum(exp(ka*(z-zm)),1),size(z,1),1)-repmat(sum(z.*exp(ka*(z-zm)),1),size(z,1),1)*ka.*exp(ka*(z-zm)))./repmat(sum(exp(ka*(z-zm)),1).^2,size(z,1),1);
end