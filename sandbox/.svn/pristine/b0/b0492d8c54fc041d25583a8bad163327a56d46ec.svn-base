function est = condEst( fileName )

  bla = load( fileName );

  M = sparse( bla(:,1)+1, bla(:,2)+1, bla(:,3) );

  est = condest( M );

  if ( nargout() == 0 ) 
    disp( est );
  endif
     

endfunction
