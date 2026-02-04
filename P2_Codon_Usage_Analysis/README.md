# P2 - Codon Usage Bias & ORF Finder

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![Biopython](https://img.shields.io/badge/Bioinformatics-Biopython-green)
![Docker](https://img.shields.io/badge/Container-Docker-2496ED)
![Status](https://img.shields.io/badge/Status-Educational-orange)

Uma ferramenta de bioinformática para análise genômica comparativa, focada em **Viés de Uso de Códons (Codon Usage Bias)** e detecção de **Open Reading Frames (ORFs)**. O pipeline automatiza o download de genomas do NCBI e gera relatórios estatísticos para Vírus, Bactérias e Eucariotos.

---

## 🎯 Objetivos de Aprendizado

Este projeto consolida os seguintes conceitos:
1.  **Código Genético:** Diferenças entre a Tabela Padrão (1) e a Tabela Bacteriana (11).
2.  **Genômica Computacional:** Algoritmos de busca de ORFs em 6 frames de leitura.
3.  **Engenharia de Dados:** Integração com APIs do NCBI (Entrez) para aquisição automática de datasets.
4.  **Estatística:** Quantificação de viés de códons (CUB) e conteúdo GC.

> **Nota Científica:** Para a discussão biológica dos resultados (Wobble Pairing, Otimização de Códons, etc.), consulte o documento [docs/final_report.md](docs/final_report.md).

---

## 🚀 Funcionalidades

- **Download Automatizado:** Baixa sequências diretamente do NCBI (ex: *E. coli*, *Lambda Phage*).
- **Detecção de ORFs:** Varre a fita *sense* e *antisense* (6 frames) e retorna coordenadas genômicas precisas.
- **Análise de Códons:** Calcula a frequência de tripletos e identifica os códons preferenciais do organismo.
- **Seleção Dinâmica de Tabela:** Aplica automaticamente a tabela de tradução correta baseada no organismo.
- **Relatórios:** Exporta JSON (dados brutos), CSV (planilha) e LaTeX (PDF).

---

## 📂 Organização do Código (MVC)

O projeto segue a arquitetura modular do repositório:

* `src/analyzer.py` **(Model):** Lógica de tradução, conversão de coordenadas AA->BP e contagem estatística.
* `src/downloader.py` **(Service):** Módulo de conexão com o NCBI Entrez.
* `src/view.py` **(View):** Visualização rica no terminal e gerador de relatórios.
* `src/main.py` **(Controller):** Orquestração do pipeline.

---

## ⚙️ Configuração

Para baixar os dados do NCBI, é necessário informar um e-mail para identificação e controle de tráfego.

1. Crie um arquivo `.env` na **raiz absoluta** do repositório:
   ```env
   ENTREZ_EMAIL=seu_email@exemplo.com
   ```
---
## 📦 Como Rodar

### Opção A: Via Docker (Recomendado)
Execute a partir da raiz do repositório para incluir as dependências compartilhadas.

```bash
# 1. Construir a imagem (a partir da raiz do repositório)
docker build -t bio-p2 -f P2_Codon_Usage_Analysis/Dockerfile .

# 2. Executar o pipeline (montando volumes para persistir dados)
docker run --rm -it \
  -v $(pwd)/P2_Codon_Usage_Analysis/data:/app/P2_Codon_Usage_Analysis/data \
  -v $(pwd)/P2_Codon_Usage_Analysis/results:/app/P2_Codon_Usage_Analysis/results \
  bio-p2

```

### Opção B: Execução Local (Python)

Utilize o ambiente virtual global na raiz.

```bash
# 1. Ative o venv (na raiz do repositório)
source venv/bin/activate

# 2. Instale as dependências específicas deste módulo
pip install -r P2_Codon_Usage_Analysis/requirements.txt

# 3. Execute o script
# (Necessário rodar da raiz para o Python encontrar o módulo 'utils')
python P2_Codon_Usage_Analysis/src/main.py

```

---

*Desenvolvido como parte do meu portfólio pessoal de estudos em Bioinformática - 2026*