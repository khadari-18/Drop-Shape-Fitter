%%% circle fitting at the apex %%%
function [r] = R_o(x,y)                                             % Andrew D. Horchler, horchler@gmail.com, Created 5-12-7
                                                                            % Revision: 1.3, 4-6-16                                                                            
%% Check inputs %%
if nargout > 4
    error('circfit:circfit:TooManyOutputs','Too many output arguments.');
end

if ~isvector(x) || ~isfloat(x) || ~isreal(x) || ~all(isfinite(x))
    error('circfit:circfit:NonFiniteRealVectorX',...
          'X must be a finite real vector of floating point numbers.');
end

if ~isvector(y) || ~isfloat(y) || ~isreal(y) || ~all(isfinite(y))
    error('circfit:circfit:NonFiniteRealVectory',...
          'Y must be a finite real vector of floating point numbers.');
end

lx=length(x);

if lx ~= length(y)
    error('circfit:circfit:LengthMismatch',...
          'The vectors X and Y must have the same length.');
end

if lx < 3
    error('circfit:circfit:Min3Points',...
          'The vectors X and Y must contain at least three points.');
end

x=x(:);
y=y(:);
% Check collinearity, assume with sufficient points, some will be non-collinear
if rank(diff([x([1:min(50,lx) 1]) y([1:min(50,lx) 1])])) == 1
    if lx <= 50 || rank(diff([x y;x(1) y(1)])) == 1
        error('circfit:circfit:Colinearity',...
             ['The points in vectors X and Y must not all be collinear, or '...
              'nearly collinear, with each other.']);
    end
end
xx=x.*x;
yy=y.*y;
xy=x.*y;
xxyy=xx+yy;
sx=sum(x);
sy=sum(y);
sxx=sum(xx);
syy=sum(yy);
sxy=sum(xy);
% Solve linear system without inverting
% a=[sx sy lx;sxy syy sy;sxx sxy sx]\[sxx+syy;sum(xxyy.*y);sum(xxyy.*x)];
[L,U]=lu([sx sy lx;sxy syy sy;sxx sxy sx]);
a=U\(L\[sxx+syy;sum(xxyy.*y);sum(xxyy.*x)]);
xc=0.5*a(1);          	% X-position of center of fitted circle
yc=0.5*a(2);          	% Y-position of center of fitted circle
r=sqrt(xc^2+yc^2+a(3));	% Radius of fitted circle

% % Set variable outputs (circle center --> xc,yc) 
% if nargout > 2
%     varargout{1}=xc;
%     varargout{2}=yc;
% end
% 
% if nargout == 2 || nargout == 4
%     varargout{nargout-1}=sqrt(mean((sqrt((x-xc).^2+(y-yc).^2)-r).^2));	% RMSE
% end

end