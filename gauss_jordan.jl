function gauss_jordan(Matriz::Matrix{T}) where {T<:Number}

  #Indices necessários para percorrer os pivores
  indiceLinha, indiceColuna = (1, 1)

  #Conversão para float para evitar exceção InexactError: Int64()
  (T <: Integer) && (Matriz = convert.(Float64, Matriz))

  # Checando se a Matriz possui uma inversa
  totalLinhas, totalColunas = size(Matriz)

  while indiceLinha<=totalLinhas && indiceColuna<=totalColunas
    
    if Matriz[indiceLinha, indiceColuna] != 0.0 
      Matriz[indiceLinha,:] = Matriz[indiceLinha,:]/Matriz[indiceLinha,indiceColuna]
      aplicaGaussJordan(Matriz, indiceLinha, indiceColuna)
    else
      trocarLinhas(indiceLinha, totalLinhas)
    end
    indiceLinha+=1
    indiceColuna+=1
  end
  return Matriz
end

function aplicaGaussJordan(Matriz, indiceLinha, indiceColuna)
  
  indicesDaColunaPivo = collect(1:length(Matriz[:,indiceColuna]))
  linhasASeremLimpas = filter((x) -> x != indiceLinha, indicesDaColunaPivo)
  for linhaASerLimpa in linhasASeremLimpas
    multiplicador = UniformScaling(Matriz[linhaASerLimpa,indiceColuna])
    Matriz[linhaASerLimpa,:] = Matriz[linhaASerLimpa,:] - (multiplicador*Matriz[indiceLinha,:])
  end
  @bp
end

function trocarLinhas(i::T, nlinha::T) where {T<:Integer}
  for n in (i+1)  #Iterando sobre as linhas acima para checar se podem ser trocadas
    if A[n,i] ≠ 0.0       # condição para trocar linha
      L = copy(A[i,:])    # copiar linha para troca
      A[i,:] = A[n,:]     # realizar troca
      A[n,:] = L
      break
    end
  end
end
