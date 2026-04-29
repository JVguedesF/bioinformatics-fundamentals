# Bioinformatics Fundamentals 🧬

![Status](https://img.shields.io/badge/Status-Study_Project-blue)
![Level](https://img.shields.io/badge/Level-Fundamentals-blue)
![Stack](https://img.shields.io/badge/Stack-Python_%7C_Biopython_%7C_Docker-2496ED)

Este repositório reúne meus **primeiros projetos de contato com a Bioinformática**. O objetivo é simples: aprender na prática os conceitos básicos de Biologia Molecular Computacional. Cada módulo é um exercício de estudo — não são ferramentas de produção.

---

## 🔬 Visão Geral dos Projetos

* **Protein & RNA Viewer (`protein-rna-viewer`):** Primeiros passos com biofísica e termodinâmica de ácidos nucleicos, predição de estruturas secundárias de RNA e visualização de estruturas proteicas.
* **Codon Bias Analyzer (`codon-bias-analyzer`):** Identificação de ORFs, tradução com tabelas genéticas e análise básica de viés de uso de códons (CUB).
* **Splicing Pattern Analyzer (`splicing-pattern-analyzer`):** Mapeamento simples de íntrons e éxons em genes eucarióticos.
* **Mobile Elements Detector (`mobile-elements-detector`):** Análise de densidade gênica, GC Skew e detecção de elementos móveis em genomas.
* **DNA Mutation Simulator (`dna-mutation-simulator`):** Simulação introdutória de replicação, sistemas de reparo e assinaturas mutacionais.

---

## 🛠️ Padrões de Engenharia

Os projetos seguem alguns princípios básicos de organização:

1. **Arquitetura MVC:** Separação entre lógica, interface e orquestração.
2. **Dockerização:** Cada módulo tem seu ambiente isolado para reprodutibilidade.
3. **Código legível:** Foco em modularidade e tipagem.

---

## 🚀 Como Executar

Cada projeto tem suas próprias instruções. Entre no diretório correspondente e consulte o `README.md` específico.

Como regra geral:
* Comandos devem ser executados a partir da raiz de cada projeto.
* Projetos que acessam o NCBI exigem um `.env` com a variável `ENTREZ_EMAIL`.

---

> [!NOTE]
> Este repositório é um portfólio de aprendizado pessoal. Os projetos foram desenvolvidos como exercícios introdutórios para exploração de conceitos básicos de Bioinformática — não se trata de software de análise científica rigorosa.
