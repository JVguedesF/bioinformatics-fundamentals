# Bioinformatics Fundamentals 🧬

![Status](https://img.shields.io/badge/Status-Active_Development-green)
![Level](https://img.shields.io/badge/Level-Fundamentals-blue)
![Stack](https://img.shields.io/badge/Stack-Python_%7C_Biopython_%7C_Docker-2496ED)

Este repositório contém a **Fase 1** da minha trilha de especialização em Bioinformática. Aqui estão implementados os conceitos fundamentais da Biologia Molecular Computacional, estruturados com padrões de engenharia de software de mercado (MVC, Docker, Clean Code).

---

## 📂 Visão Geral dos Projetos

| ID | Projeto | Foco Biológico | Stack Principal |
| :--- | :--- | :--- | :--- |
| **[P1](./P1_Structural_Profiling_Viz)** | **Structural Profiling** | Biofísica, Termodinâmica (DNA/RNA) e Estrutura 3D | `Biopython` `ViennaRNA` `PyMOL` |
| **[P2](./P2_Codon_Usage_Analysis)** | **Codon Usage Analysis** | Genômica Comparativa, ORFs e Viés de Códons | `NCBI Entrez` `Rich` `Biopython` |
| **[P3](./P3_Eukaryotic_Splicing_Dynamics)** | **Eukaryotic Splicing** | Transcriptômica e Isoformas | *(Em breve)* |
| **[P4](./P4_Mobile_Elements_Genomics)** | **Mobile Elements** | Genomas Organelares e Transposons | *(Em breve)* |
| **[P5](./P5_Replication_Mutation_Sim)** | **Replication Sim** | Dinâmica de Replicação e Reparo | *(Em breve)* |

---

## 🛠️ Padrões de Engenharia

Para garantir reprodutibilidade e organização, todos os projetos seguem uma arquitetura padronizada:

1.  **Monorepo Modular:** Código compartilhado reside na pasta `utils/`, evitando duplicação.
2.  **Arquitetura MVC:** Separação estrita entre Lógica Biológica (`Model`), Interface CLI (`View`) e Orquestração (`Controller`).
3.  **Docker First:** Cada projeto possui seu próprio container isolado para resolver dependências de sistema.

---

## 🚀 Como Executar

Cada projeto funciona como um módulo independente com sua própria documentação e container Docker.

1.  **Escolha o projeto** na tabela acima.
2.  Acesse a pasta correspondente (ex: `cd P1_Structural_Profiling_Viz`).
3.  Siga as instruções do **`README.md` local** para construir a imagem Docker ou rodar o script Python.

---

## ⚙️ Configuração

Para baixar os dados do NCBI, é necessário informar um e-mail para identificação e controle de tráfego.

1. Crie um arquivo `.env` na **raiz absoluta** do repositório:
   ```env
   ENTREZ_EMAIL=seu_email@exemplo.com

---

## 📦 Instalação Local (Desenvolvimento)

Se preferir rodar sem Docker (via IDE), configure o ambiente virtual na raiz para que todos os projetos compartilhem as dependências base:

```bash
# 1. Criar e ativar venv na raiz
python3 -m venv venv
source venv/bin/activate  # ou venv\Scripts\activate no Windows

# 2. Instalar dependências dos módulos
pip install -r P1_Structural_Profiling_Viz/requirements.txt
pip install -r P2_Codon_Usage_Analysis/requirements.txt

# 3. Rodar um projeto (sempre a partir da raiz)
# Exemplo: Rodando o P2
python P2_Codon_Usage_Analysis/src/main.py