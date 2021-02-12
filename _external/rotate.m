function out = rotate(img, angle)

[rowsi,colsi,z]= size(img);


switch mod(angle, 360)
    % Special cases
    case 0
        out = img;
    case 90
        out = rot90(img);
    case 180
        out = img(end:-1:1, end:-1:1);
    case 270
        out = rot90(img(end:-1:1, end:-1:1));
        
        % General rotations
    otherwise
        rads=2*pi*angle/360;
        
        %calculating array dimesions such that  rotated image gets fit in it exactly.
        % we are using absolute so that we get  positve value in any case ie.,any quadrant.
        
        rowsf=ceil(rowsi*abs(cos(rads))+colsi*abs(sin(rads)));
        colsf=ceil(rowsi*abs(sin(rads))+colsi*abs(cos(rads)));
        
        % define an array withcalculated dimensionsand fill the array  with zeros ie.,black
        C = zeros([rowsf colsf ]);
        
        %calculating center of original and final image
        xo=ceil(rowsi/2);
        yo=ceil(colsi/2);
        
        midx=ceil((size(C,1))/2);
        midy=ceil((size(C,2))/2);
        
        % in this loop we calculate corresponding coordinates of pixel of A
        % for each pixel of C, and its intensity will be  assigned after checking
        % weather it lie in the bound of A (original image)
        for i=1:size(C,1)
            for j=1:size(C,2)
                
                x= (i-midx)*cos(rads)+(j-midy)*sin(rads);
                y= -(i-midx)*sin(rads)+(j-midy)*cos(rads);
                x=round(x)+xo;
                y=round(y)+yo;
                
                if (x>=1 && y>=1 && x<=size(img,1) &&  y<=size(img,2) )
                    C(i,j,:)=img(x,y,:);
                end
                
            end
        end
        
        out = tools.crop(C,size(img));
end