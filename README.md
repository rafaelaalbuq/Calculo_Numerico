# Simulador de Zeros de Funções — Cálculo Numérico (UNIVASF 2026.2)

[![Python 3](https://img.shields.io/badge/Python-3.11%2B-blue.svg)](https://www.python.org/)
[![Status](https://img.shields.io/badge/Status-Concluído-brightgreen.svg)]()
[![Disciplina](https://img.shields.io/badge/UNIVASF-Cálculo%20Numérico-navy.svg)]()

> **Projeto Prático da Unidade 1 — Atividade 2: Comparativo de Métodos para Resolução de Equações Algébricas e Transcendentes**  
> Fundação Universidade Federal do Vale do São Francisco (UNIVASF) — Semestre 2026.2

---

## 👥 Equipe
- **Marcos Anderson P. de Carvalho**
- **Pedro de Souza Torres**
- **Rafaela Nunes Albuquerque**

---

## 📌 Visão Geral do Projeto

Este projeto consiste em um simulador numérico completo desenvolvido em **Python 3** para determinação e comparação de raízes (zeros) de equações algébricas e transcendentes. A aplicação implementa quatro dos métodos iterativos clássicos do Cálculo Numérico:

1. **Método da Bissecção** *(Método de Intervalo / Fechado)*
2. **Método da Falsa Posição** *(Método de Intervalo / Fechado)*
3. **Método de Newton-Raphson** *(Método Aberto)*
4. **Método da Secante** *(Método Aberto)*

Além do cálculo das raízes, o software realiza:
- **Medição precisa de tempo** de execução via `time.perf_counter` (média estabilizada de 500 execuções, medida em microssegundos $\mu s$);
- **Contabilização rigorosa de iterações** até atingir o critério de parada ($|f(x_k)| < \varepsilon$);
- **Tratamento preventivo de exceções e patologias numéricas** (divisão por zero, derivada nula, ausência de mudança de sinal e extrapolação de domínio físico);
- **Interface flexível** que permite ao usuário inserir qualquer nova equação por texto;
- **Plotagem e exportação de gráficos** comparativos (função + zeros + curva de convergência em escala logarítmica).

---

## 🚀 Como Executar

### 1. Clonar o Repositório
```bash
git clone https://github.com/rafaelaalbuq/Calculo_Numerico.git
cd Calculo_Numerico
```

### 2. Instalar as Dependências
Recomenda-se o uso de um ambiente virtual ou instalação direta das bibliotecas listadas no `requirements.txt`:
```bash
pip install -r requirements.txt
```
*(Dependências: `matplotlib`, `numpy`, `reportlab`)*

### 3. Modos de Execução

#### A) Modo Interativo (Menu no Terminal)
Permite navegar entre os casos cadastrados, testar funções customizadas e realizar novos cálculos:
```bash
python simulador_zeros.py
```

#### B) Modo Automatizado em Lote (Gera Gráficos)
Roda todos os casos de teste e problemas de engenharia e salva as imagens PNG:
```bash
python simulador_zeros.py --auto
```

#### C) Modo Sem Renderização Gráfica
Ideal para execução rápida em servidores ou terminais sem interface gráfica:
```bash
python simulador_zeros.py --auto --sem-grafico
```

#### D) Gerar o Relatório Técnico em PDF
Gera o relatório oficial de 4 páginas com todas as análises, tabelas e gráficos incorporados:
```bash
python gerar_relatorio.py
```

---

## 📊 Casos Analisados e Resultados Obtidos

Todos os testes foram executados com tolerância $\varepsilon = 10^{-6}$ e limite máximo de 100 iterações.

### 1. Exemplo de Demonstração: $f(x) = \cos(x) - x = 0$
- **Intervalo:** $[0, 1]$ | **Newton:** $x_0 = 0.5$ | **Secante:** $x_0 = 0, x_1 = 1$
- **Raiz Analítica de Referência:** $x^* \approx 0.7390851332$ (Número de Dottie)

| Método | Raiz Encontrada | Iterações | Tempo Médio ($\mu s$) | $|f(\text{raiz})|$ | Erro vs Ref. | Status |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Bissecção** | 0.7390851974 | 20 | 7.69 | $1.075 \times 10^{-7}$ | $6.42 \times 10^{-8}$ | OK |
| **Falsa Posição** | 0.7390847824 | 5 | 4.41 | $5.870 \times 10^{-7}$ | $3.51 \times 10^{-7}$ | OK |
| **Newton-Raphson** | 0.7390851339 | 3 | 2.91 | $1.181 \times 10^{-9}$ | $7.21 \times 10^{-10}$ | OK |
| **Secante** | 0.7390851121 | 4 | 2.70 | $3.529 \times 10^{-8}$ | $2.11 \times 10^{-8}$ | OK |

---

### 2. Problema 1: Resfriamento de um Componente Eletrônico
- **Equação:** $f(t) = 20 + 90e^{-0.15t} - 50 = 0$, com $t \ge 0$
- **Derivada:** $f'(t) = -13.5e^{-0.15t}$
- **Solução Exata:** $t = \frac{\ln(3)}{0.15} \approx 7.3240819245\text{ minutos}$ (~7 min 19,4 s)
- **Intervalo:** $[0, 10]$ min | **Newton:** $t_0 = 0$ | **Secante:** $t_0 = 0, t_1 = 10$

| Método | Raiz (min) | Iterações | Tempo Médio ($\mu s$) | $|f(\text{raiz})|$ | Erro vs Ref. | Status |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Bissecção** | 7.3240819573 | 25 | 12.37 | $1.480 \times 10^{-7}$ | $3.29 \times 10^{-8}$ | OK |
| **Falsa Posição** | 7.3240820805 | 21 | 12.00 | $7.023 \times 10^{-7}$ | $1.56 \times 10^{-7}$ | OK |
| **Newton-Raphson** | 7.3240819244 | 5 | 5.72 | $3.966 \times 10^{-10}$ | $8.81 \times 10^{-11}$ | OK |
| **Secante** | 7.3240819245 | 6 | 5.13 | $2.134 \times 10^{-10}$ | $4.74 \times 10^{-11}$ | OK |

---

### 3. Problema 2: Carregamento de um Capacitor em Circuito RC
- **Equação:** $f(t) = 5(1 - e^{-t}) - 3.8 = 0$, com $t \ge 0$
- **Derivada:** $f'(t) = 5e^{-t}$
- **Solução Exata:** $t = -\ln(0.24) \approx 1.4271163556\text{ segundos}$

#### Caso A: Intervalo [0, 5] (Chutes nas Extremidades)
- A **Secante** com $t_0=0$ e $t_1=5$ gera na 2ª iteração $t_2 \approx -13.18\text{ s} < 0$, violando a restrição física de causalidade ($t \ge 0$). O simulador intercepta a falha por extrapolação de domínio e encerra com segurança antes que ocorra overflow numérico.

#### Caso B: Intervalo Alternativo Refinado [0, 3]
| Método | Raiz (s) | Iterações | Tempo Médio ($\mu s$) | $|f(\text{raiz})|$ | Erro vs Ref. | Status |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Bissecção** | 1.4271168709 | 21 | 9.89 | $6.183 \times 10^{-7}$ | $5.15 \times 10^{-7}$ | OK |
| **Falsa Posição** | 1.4271170326 | 25 | 13.52 | $8.124 \times 10^{-7}$ | $6.77 \times 10^{-7}$ | OK |
| **Newton-Raphson** | 1.4271163488 | 5 | 6.07 | $8.179 \times 10^{-9}$ | $6.82 \times 10^{-9}$ | OK |
| **Secante** | 1.4271163561 | 9 | 6.15 | $5.550 \times 10^{-10}$ | $4.63 \times 10^{-10}$ | OK |

---

## 🛡️ Tratamento de Falhas e Exceções Numéricas

O simulador detecta preventivamente as principais anomalias do Cálculo Numérico:
1. **Ausência de Troca de Sinal:** Verifica se $f(a) \cdot f(b) > 0$ antes de iniciar a iteração dos métodos fechados;
2. **Derivada Nula ($f'(x_k) = 0$):** Intercepta tangentes horizontais no método de Newton-Raphson evitando erro de divisão por zero;
3. **Divisão por Zero na Secante ($f(x_k) = f(x_{k-1})$):** Detecta secantes horizontais e reporta o passo da falha;
4. **Divergência / Ciclo Infinito:** Impede travamento por meio do teto configurável `max_iter`;
5. **Extrapolação de Domínio Físico:** Garante que passos de métodos abertos não transgridam intervalos com restrições físicas (ex.: $t \ge 0$).

---

## 📁 Estrutura de Arquivos

```text
Cal_num/
├── simulador_zeros.py       # Código-fonte principal do simulador
├── gerar_relatorio.py       # Script gerador do relatório técnico em PDF
├── relatorio_zeros_univasf.pdf # Relatório final completo formatado
├── requerimentos.pdf        # Documento oficial com as diretrizes do projeto
├── requirements.txt         # Lista de dependências Python
├── README.md                # Documentação técnica do repositório
├── .gitignore               # Arquivos ignorados pelo controle de versão
└── grafico_*.png            # Gráficos de funções e curvas de convergência
```

---

## 📄 Licença e Uso Acadêmico
Projeto desenvolvido exclusivamente para fins didáticos e avaliativos na disciplina de Cálculo Numérico da **Universidade Federal do Vale do São Francisco (UNIVASF)**.

