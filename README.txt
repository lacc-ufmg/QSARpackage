1. Criar o ambiente virtual

conda env create -f environment.yml

2. Para ativar o ambiente virtual criado quando abrir um terminal

conda activate QSAR

3. Para desativar o ambiente virtual

conda deactivate

4. Para rodar o programa

4.1. De dentro do diretório QSARpackage

4.1.2. Usando arquivo contendo alinhamentos

python runMD.py -m diretorio -e extensão -a alinhamentos -p sondas -d step

-diretorio: diretório contendo arquivos com as moléculas. Pode ser qualquer formato de arquivo que o openbabel entenda

-extensão: extensão dos arquivo com as moléculas (mol2, por exemplo)

-alinhamentos: arquivo csv, com a primeira coluna contendo os nomes das moléculas e a segunda contendo os átomos para os alinhamentos. As colunas devem ser separadas por ;

-sondas: sondas utilizadas para calcular os descritores com o LQTAgrid (NH3+, por exemplo)

-step: Tamanho do passo que a sonda percorre em angstrom.

Exemplo:

python runMD.py -m ~/teste -e mol2 -a alinhamentos.csv -p NH3+ -d 1

4.1.3. Usando SMARTS para alinhar

python runMD.py -m diretorio -e extensão -s smarts -p sondas -d step

-smarts: código SMARTS comum aos átomos que serão usados para o alinhamento.

Exemplo:

python runMD.py -m ~/teste -e mol2 -s *~2~*~*~1~*(~*~*~*~*~1)~*~2 -p NH3+ -d 1

4.2. Rodando dentro do diretório onde estão localizadas as moléculas

Exemplo:

python caminho_para_pasta_QSARpackage/QSARpackage/runMD.py -m . -e mol2 -a alinhamentos.csv -p NH3+ -d 1

ou

python caminho_para_pasta_QSARpackage/QSARpackage/runMD.py -m . -e mol2 -s *~2~*~*~1~*(~*~*~*~*~1)~*~2 -p NH3+ -d 1

As opções relativas aos alinhamentos e LQTAGrid não são  obrigatórias. O alnhamento pode ser rodado posteriormente com o programa runAlignment.py e o LQTAGrid pode ser rodado com o programa runLQTAGrid.py

5. Para gerar apenas os descritores 3D ou 4D

5.1. De dentro do diretório QSARpackage python runLQTAGrid.py -m diretorio -e extensão -a NH3+ -s step -o diretorio/matriz

-diretorio: diretório contendo arquivos com as moléculas. Pode ser qualquer formato de arquivo que o openbabel entenda

-extensão: extensão dos arquivo com as moléculas (mol2, por exemplo)

-alinhamentos: arquivo csv, com a primeira coluna contendo os nomes das moléculas e a segunda contendo os átomos para os alinhamentos. As colunas devem ser separadas por ;

-sondas: sondas utilizadas para calcular os descritores com o LQTAgrid (NH3+, por exemplo)

-step: Tamanho do passo que a sonda percorre em angstrom.


