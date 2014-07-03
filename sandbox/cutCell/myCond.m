function c = myCond( name )

bla = load( name );

M = sparse( bla(:,1)+1, bla(:,2)+1, bla(:,3) );

%%c = cond( M );
c = condest(M);
display(c)

endfunction


%% -- Call silently from the shell with:
%%
%% echo -n `octave --silent --eval "myCond('matrix');"`; echo -e ''
%%
