tic;
% n �������л��ֵĸ�������������[0,1]*[0,1]��n^2������,canΪ�����еĲ���k
%����������1,2,3����ʱ������ʾһ�������εĸ�������
% s��һ��n^2*10�Ĺ�������s(i,1)��ʾ��i�������Σ�s(i,2),s(i,3),s(i,4)�ֱ��ʾ��i�������ε�1,2,3����Ӧ�Ķ���
% s(i,5),s(i,6)��s(i,7)��s(i,8)��s(i,9)��s(i,10)�ֱ��ʾ����1,2,3�����������
%���ɹ�������s
%A���ܸվ���
%�������ű���x y
n=20;
can=20;
s=zeros(n^2,10);
h=1/n;
st=1/(2*n^2);
A=zeros((n+1)^2,(n+1)^2);
syms x y;
for  k=1:1:2*n^2
        s(k,1)=k;
        q=fix(k/(2*n));
        r=mod(k,(2*n));
        if (r~=0)
            r=r;
        else r=2*n;q=q-1;
        end
        if (r<=n)
            s(k,2)=q*(n+1)+r;
            s(k,3)=q*(n+1)+r+1;
            s(k,4)=(q+1)*(n+1)+r+1;
            s(k,5)=(r-1)*h;
            s(k,6)=q*h;
            s(k,7)=r*h;
            s(k,8)=q*h;
            s(k,9)=r*h;
            s(k,10)=(q+1)*h;
        else
            s(k,2)=q*(n+1)+r-n;
            s(k,3)=(q+1)*(n+1)+r-n+1;
            s(k,4)=(q+1)*(n+1)+r-n;
            s(k,5)=(r-n-1)*h;
            s(k,6)=q*h;
            s(k,7)=(r-n)*h;
            s(k,8)=(q+1)*h;
            s(k,9)=(r-n-1)*h;
            s(k,10)=(q+1)*h;
        end
end
%�������ɻ�����L(i)��ʾ��i���㶥��Ļ�����
%���ɵ��վ���d
%���ɵ��վ��󲢽�������ܸپ���
d=zeros(3,3);
B=zeros((n+1)^2,1);
b=zeros(3,1);
%����A���ܸ�
for   k=1:1:2*n^2
     L(1)=(1/(2*st))*((s(k,7)*s(k,10)-s(k,9)*s(k,8))+(s(k,8)-s(k,10))*x+(s(k,9)-s(k,7))*y);
     L(2)=(1/(2*st))*((s(k,9)*s(k,6)-s(k,5)*s(k,10))+(s(k,10)-s(k,6))*x+(s(k,5)-s(k,9))*y);
     L(3)=(1/(2*st))*((s(k,5)*s(k,8)-s(k,7)*s(k,6))+(s(k,6)-s(k,8))*x+(s(k,7)-s(k,5))*y);
     for i=1:1:3
        for j=i:3
            d(i,j)=int(int(((((diff(L(i),x))*(diff(L(j),x)))+((diff(L(i),y))*(diff(L(j),y))))-((can^2)*L(i)*L(j))),x,0,1),y,0,1);
            d(j,i)=d(i,j);
        end
     end   
     for i=1:1:3
         for j=1:1:3
             A(s(k,(i+1)),s(k,(j+1)))=A(s(k,(i+1)),s(k,(j+1)))+d(i,j);
         end
     end
     for i=1:1:3
         b(i)=int(int((L(i)),x,0,1),y,0,1);
         B(s(k,(i+1)),1)=B(s(k,(i+1)),1)+b(i);
     end
     
end
%������ݱ߽��������������Ԫ������Mx=B,��α߽�����Լ���˺ܶ���
M=zeros((n+1)^2,n^2);
j=n^2;
for i=(n^2+n):-1:1
    if ((mod(i,(n+1)))~=1)
        M(:,j)=A(:,i);
        j=j-1;
    else continue
    end
end
%preanswer��δ֪���ֵ�ǣ�n+1��^2*(n^2)��
preanswer=M\B;
%�õ����нڵ��ֵ
answer=zeros((n+1)^2,1);
j=1;
for i=1:1:(n^2+n)
    if ((mod(i,(n+1)))~=1)
        answer(i)=preanswer(j);
        j=j+1;
    else answer(i)=0;
    end
end
%���������ͼ
Z=zeros((n+1),(n+1));
for i=1:1:(n+1)^2
      s=fix(i/(n+1))+1;
      r=mod(i,(n+1));
      if(r==0)
          r=n+1;
          s=s-1;
      else
      end
      Z(r,s)=answer(i);
end
[X,Y]=meshgrid(1:-h:0,0:h:1);
surf(X,Y,Z);



toc;
t=toc;

     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
      