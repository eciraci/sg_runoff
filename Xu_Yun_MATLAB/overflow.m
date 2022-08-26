function [ixoverflow,iyoverflow]=overflow(ix,iy)
global phi isold

for ixn=ix-1:ix+1
  for iyn=iy-1:iy+1
    if  ~isold(ixn,iyn)
      isold(ixn,iyn)=logical(1);
%----
if phi(ixn,iyn)<phi(ix,iy)
  ixoverflow=ixn;iyoverflow=iyn;
elseif phi(ixn,iyn) == phi(ix,iy) 
    [ixoverflow,iyoverflow] = overflow(ixn,iyn) ;
end
%----
    end
  end
end
