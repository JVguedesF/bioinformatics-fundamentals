# P1 - Structural Profiling & Molecular Visualization

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![Biopython](https://img.shields.io/badge/Bioinformatics-Biopython-green)
![Docker](https://img.shields.io/badge/Container-Docker-2496ED)
![Status](https://img.shields.io/badge/Status-Educational-orange)

Projeto introdutório de bioinformática focado na manipulação básica de sequências biológicas. O objetivo deste módulo é automatizar cálculos fundamentais para DNA, RNA e Proteínas, servindo como base para projetos mais complexos e exercitando boas práticas de organização de código em Python.

---

## 🎯 Objetivos de Aprendizado (Tech & Bioinfo)

Este projeto foi desenvolvido para consolidar:
1. **Domínio do Biopython:** Manipulação de arquivos `.fasta` e cálculos físico-químicos (Peso Molecular, GC%, pI).
2. **Integração de Software Externo:** Uso de wrappers para rodar o *ViennaRNA* (para estrutura de RNA) dentro do Python.
3. **Padrão MVC:** Separação da lógica de cálculo (Model) da interface de terminal (View).
4. **Conteinerização:** Uso de Docker para evitar problemas de dependência de bibliotecas C (como o RNAfold).

> **Nota Científica:** Para a discussão biológica dos resultados, validação dos dados e gráficos de estrutura 3D, consulte o documento [docs/final_report.md](docs/final_report.md).

---

## ⚙️ Funcionalidades

O script detecta automaticamente o tipo de sequência no arquivo FASTA e aplica o pipeline adequado:

* **DNA Pipeline:** * Cálculo de GC Content (%).
  * Estimativa de Temperatura de Melting (Tm).
* **RNA Pipeline:** * Predição de Estrutura Secundária e Energia Livre Mínima (MFE) via *ViennaRNA*.
* **Protein Pipeline:** * Cálculo de Peso Molecular (Da) e Ponto Isoelétrico (pI).
  * Índice de Instabilidade e GRAVY (hidrofobicidade).

**Saída:** Visualização interativa no terminal (via `rich`) e exportação para arquivos JSON e CSV na pasta `results/`.

---

## 📂 Organização do Código (MVC)

O código foi refatorado para evitar scripts monolíticos, dividindo responsabilidades:

* `src/analyzer.py` **(Model):** Contém toda a lógica biológica e chamadas ao Biopython. Não possui funções de `print`.
* `src/view.py` **(View):** Responsável exclusivamente por renderizar tabelas e textos coloridos no terminal usando a biblioteca `rich`.
* `src/main.py` **(Controller):** Orquestra o fluxo. Lê os arquivos `data/sequences/`, envia para o Analyzer e manda os resultados para a View.

---

## 🚀 Como Rodar

### Opção A: Usando Docker (Recomendado)
O Docker garante que as dependências do ViennaRNA funcionem independentemente do seu sistema operacional.

```bash
# 1. Construir a imagem (a partir da raiz do repositório)
docker build -t bioinfo-p1 -f P1_Structural_Profiling_Viz/Dockerfile .

# 2. Executar o container montando as pastas locais
docker run --rm \
  -v $(pwd)/P1_Structural_Profiling_Viz/data:/app/P1_Structural_Profiling_Viz/data \
  -v $(pwd)/P1_Structural_Profiling_Viz/results:/app/P1_Structural_Profiling_Viz/results \
  bioinfo-p1

```

### Opção B: Execução Local (Python)
Este projeto utiliza o ambiente virtual configurado na raiz do repositório.

```bash
# 1. Ative o ambiente virtual (a partir da raiz do repositório)
source venv/bin/activate

# 2. Entre na pasta do projeto P1
cd P1_Structural_Profiling_Viz

# 3. Instale as dependências (caso ainda não tenha feito)
pip install -r requirements.txt

# 4. Execute o script
python src/main.py
```

---

## 🛠️ Stack Tecnológica

* **Linguagem:** Python 3.10
* **Bibliotecas Biológicas:** Biopython, ViennaRNA (RNAfold)
* **CLI/Visualização:** Rich

---

*Desenvolvido como parte do meu portfólio pessoal de estudos em Bioinformática - 2026*