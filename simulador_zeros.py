#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulador de Resolução de Equações Algébricas e Transcendentes
Disciplina: Cálculo Numérico - UNIVASF - Projeto Unidade 1 (Zeros)

Métodos implementados:
    1. Bissecção          (método de intervalo / fechado)
    2. Falsa Posição      (método de intervalo / fechado)
    3. Newton-Raphson     (método aberto)
    4. Secante            (método aberto)

Uso:
    python simulador_zeros.py          -> menu interativo
    python simulador_zeros.py --auto   -> roda todos os casos exigidos e salva os gráficos
    python simulador_zeros.py --auto --sem-grafico
"""

import math
import sys
import time
from dataclasses import dataclass, field

# Garante compatibilidade de saída UTF-8 no terminal Windows (evita UnicodeEncodeError em cp1252)
if hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


# Parâmetros globais (mesmos para todos os métodos => comparação justa)

EPS_PADRAO = 1e-6        # tolerância: o método para quando |f(x)| fica abaixo deste valor
MAX_ITER_PADRAO = 100    # limite de iterações para evitar que um método rode indefinidamente
REPETICOES_TEMPO = 500   # repetições usadas para calcular uma média de tempo mais estável



# Estrutura que guarda o resultado de cada método

@dataclass
class Resultado:
    """Reúne a resposta e as informações de execução de um método numérico."""
    metodo: str
    raiz: float = float("nan")
    iteracoes: int = 0
    tempo_us: float = float("nan")      # tempo médio em microssegundos
    residuo: float = float("nan")       # |f(raiz)|
    convergiu: bool = False
    mensagem: str = ""                  # explica falhas ou o critério atingido
    historico: list = field(default_factory=list)  # |f(x_k)| a cada iteração



# MÉTODOS NUMÉRICOS
# Critério de parada (todos): |f(x_k)| < eps  

def bisseccao(f, a, b, eps=EPS_PADRAO, max_iter=MAX_ITER_PADRAO):
    """Procura uma raiz em [a, b], mantendo sempre um intervalo com troca de sinal.

    A cada passo usa o ponto médio; a raiz fica na metade em que f muda de sinal.
    """
    r = Resultado("Bissecção")
    fa, fb = f(a), f(b)
    # Sem troca de sinal, este método não tem como garantir que há uma raiz no intervalo.
    if fa * fb > 0:
        r.mensagem = "FALHA: f(a) e f(b) têm o mesmo sinal (sem mudança de sinal no intervalo)."
        return r
    x_ant = a
    for k in range(1, max_iter + 1):
        # Testa o ponto médio e guarda o resíduo para o gráfico de convergência.
        m = (a + b) / 2.0
        fm = f(m)
        r.historico.append(abs(fm))
        if abs(fm) < eps:
            r.raiz, r.iteracoes, r.residuo, r.convergiu = m, k, abs(fm), True
            r.mensagem = "Convergiu."
            return r
        if fa * fm < 0:      # raiz está em [a, m]
            b, fb = m, fm
        else:                # raiz está em [m, b]
            a, fa = m, fm
        x_ant = m
    r.raiz, r.iteracoes, r.residuo = m, max_iter, abs(fm)
    r.mensagem = "FALHA: número máximo de iterações atingido."
    return r


def falsa_posicao(f, a, b, eps=EPS_PADRAO, max_iter=MAX_ITER_PADRAO):
    """Estima a raiz pela interseção da reta entre os extremos com o eixo x.

    Como na bissecção, conserva o subintervalo em que f troca de sinal.
    """
    r = Resultado("Falsa Posição")
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        r.mensagem = "FALHA: f(a) e f(b) têm o mesmo sinal (sem mudança de sinal no intervalo)."
        return r
    x_ant = a
    for k in range(1, max_iter + 1):
        # Esta fórmula calcula onde a reta que liga (a, f(a)) e (b, f(b)) cruza y=0.
        if fb - fa == 0:
            r.mensagem = "FALHA: divisão por zero (f(b) - f(a) = 0)."
            return r
        x = (a * fb - b * fa) / (fb - fa)
        fx = f(x)
        r.historico.append(abs(fx))
        if abs(fx) < eps:
            r.raiz, r.iteracoes, r.residuo, r.convergiu = x, k, abs(fx), True
            r.mensagem = "Convergiu."
            return r
        if fa * fx < 0:
            b, fb = x, fx
        else:
            a, fa = x, fx
        x_ant = x
    r.raiz, r.iteracoes, r.residuo = x, max_iter, abs(fx)
    r.mensagem = "FALHA: número máximo de iterações atingido."
    return r


def newton_raphson(f, df, x0, eps=EPS_PADRAO, max_iter=MAX_ITER_PADRAO, dominio=None):
    """Aplica Newton-Raphson usando f e sua derivada df a partir de x0.

    O argumento opcional dominio limita os valores aceitos para os novos iterados.
    """
    r = Resultado("Newton-Raphson")
    x = x0
    try:
        for k in range(1, max_iter + 1):
            fx, dfx = f(x), df(x)
            # A fórmula exige dividir por f'(x); se a derivada for zero, não há passo válido.
            if dfx == 0:
                r.mensagem = f"FALHA: derivada nula em x = {x:.6g} (iteração {k})."
                return r
            x_novo = x - fx / dfx  # interseção da tangente em x com o eixo horizontal
            if dominio and not (dominio[0] <= x_novo <= dominio[1]):
                r.raiz, r.iteracoes = x_novo, k
                r.mensagem = (f"FALHA: extrapolação do domínio (x_{k} = {x_novo:.6g} "
                              f"fora de [{dominio[0]}, {dominio[1]}]).")
                return r
            f_novo = f(x_novo)
            r.historico.append(abs(f_novo))
            if abs(f_novo) < eps:
                r.raiz, r.iteracoes, r.residuo, r.convergiu = x_novo, k, abs(f_novo), True
                r.mensagem = "Convergiu."
                return r
            x = x_novo
    except (OverflowError, ZeroDivisionError) as e:
        r.mensagem = f"FALHA: erro numérico ({e})."
        return r
    r.raiz, r.iteracoes, r.residuo = x, max_iter, abs(f(x))
    r.mensagem = "FALHA: número máximo de iterações atingido."
    return r


def secante(f, x0, x1, eps=EPS_PADRAO, max_iter=MAX_ITER_PADRAO, dominio=None):
    """Aproxima a raiz usando dois pontos iniciais, sem precisar calcular f' explicitamente."""
    r = Resultado("Secante")
    f0, f1 = f(x0), f(x1)
    try:
        for k in range(1, max_iter + 1):
            # Se os valores de f forem iguais, a reta secante é horizontal e não cruza o eixo x.
            if f1 - f0 == 0:
                r.mensagem = f"FALHA: divisão por zero (f(x_k) = f(x_(k-1))) na iteração {k}."
                return r
            x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
            if dominio and not (dominio[0] <= x2 <= dominio[1]):
                r.raiz, r.iteracoes = x2, k
                r.mensagem = (f"FALHA: extrapolação do domínio (x = {x2:.6g} "
                              f"fora de [{dominio[0]}, {dominio[1]}]).")
                return r
            f2 = f(x2)
            r.historico.append(abs(f2))
            if abs(f2) < eps:
                r.raiz, r.iteracoes, r.residuo, r.convergiu = x2, k, abs(f2), True
                r.mensagem = "Convergiu."
                return r
            x0, f0, x1, f1 = x1, f1, x2, f2
    except (OverflowError, ZeroDivisionError) as e:
        r.mensagem = f"FALHA: erro numérico ({e})."
        return r
    r.raiz, r.iteracoes, r.residuo = x1, max_iter, abs(f1)
    r.mensagem = "FALHA: número máximo de iterações atingido."
    return r



# Medição de tempo e execução comparativa

def medir(execucao, repeticoes=REPETICOES_TEMPO):
    """Executa o método uma vez (resultado) e depois `repeticoes` vezes para obter
    o tempo médio em microssegundos (uma única execução é rápida demais para medir)."""
    resultado = execucao()  # guarda uma execução real para exibir raiz, iterações e histórico
    t0 = time.perf_counter()
    for _ in range(repeticoes):
        execucao()
    resultado.tempo_us = (time.perf_counter() - t0) / repeticoes * 1e6
    return resultado


def comparar(f, df, a, b, x0_newton, x0_sec, x1_sec, eps=EPS_PADRAO,
             max_iter=MAX_ITER_PADRAO, dominio=None):
    """Executa os quatro métodos com os mesmos limites e devolve seus resultados."""
    return [
        medir(lambda: bisseccao(f, a, b, eps, max_iter)),
        medir(lambda: falsa_posicao(f, a, b, eps, max_iter)),
        medir(lambda: newton_raphson(f, df, x0_newton, eps, max_iter, dominio)),
        medir(lambda: secante(f, x0_sec, x1_sec, eps, max_iter, dominio)),
    ]


def imprimir_tabela(titulo, resultados, raiz_ref=None):
    """Mostra raiz, custo, resíduo e estado de cada método em formato de tabela."""
    print("\n" + "=" * 92)
    print(titulo)
    print("=" * 92)
    print(f"{'Método':<16}{'Raiz':>16}{'Iter.':>7}{'Tempo (µs)':>13}{'|f(raiz)|':>14}{'Erro vs ref.':>14}  Status")
    print("-" * 92)
    for r in resultados:
        if math.isnan(r.raiz) or not r.convergiu:
            print(f"{r.metodo:<16}{'---':>16}{r.iteracoes:>7}{'---':>13}{'---':>14}{'---':>14}  {r.mensagem}")
            continue
        erro = f"{abs(r.raiz - raiz_ref):.2e}" if raiz_ref is not None else "---"
        print(f"{r.metodo:<16}{r.raiz:>16.10f}{r.iteracoes:>7}{r.tempo_us:>13.2f}"
              f"{r.residuo:>14.3e}{erro:>14}  OK")
    print("-" * 92)
    print(f"ε = {EPS_PADRAO:g} | máx. iterações = {MAX_ITER_PADRAO} | "
          f"tempo = média de {REPETICOES_TEMPO} execuções")



# Gráficos (bônus): função + zeros e convergência de |f(x_k)|

def plotar(titulo, f, a, b, resultados, arquivo, mostrar=True):
    """Desenha a função e a convergência; salva o PNG e, se pedido, abre a janela."""
    import matplotlib
    if not mostrar:
        # Backend sem janela: permite salvar gráficos em execuções automáticas/headless.
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    xs = np.linspace(a, b, 600)  # pontos usados para desenhar a curva da função
    ax1.plot(xs, [f(x) for x in xs], "b-", label="f(x)")
    ax1.axhline(0, color="k", lw=0.8)
    marcadores = {"Bissecção": "o", "Falsa Posição": "s", "Newton-Raphson": "^", "Secante": "D"}
    for i, r in enumerate(resultados):
        if r.convergiu:
            ax1.plot(r.raiz, 0, marcadores[r.metodo], ms=11 - 2 * i, mfc="none", mew=2,
                     label=f"{r.metodo}: {r.raiz:.6f}")
    ax1.set_title(f"{titulo}\nFunção e raízes encontradas")
    ax1.set_xlabel("x"); ax1.set_ylabel("f(x)"); ax1.grid(alpha=0.3); ax1.legend(fontsize=8)

    for r in resultados:
        if r.historico:
            ax2.semilogy(range(1, len(r.historico) + 1),
                         [max(h, 1e-17) for h in r.historico],
                         marker=marcadores[r.metodo], label=r.metodo)
    ax2.axhline(EPS_PADRAO, color="r", ls="--", label=f"ε = {EPS_PADRAO:g}")
    ax2.set_title("Convergência: |f(x_k)| por iteração")
    ax2.set_xlabel("Iteração k"); ax2.set_ylabel("|f(x_k)| (escala log)")
    ax2.grid(alpha=0.3, which="both"); ax2.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(arquivo, dpi=150)  # salva no diretório atual ou no caminho informado
    print(f"[gráfico salvo em: {arquivo}]")
    if mostrar:
        plt.show()
    plt.close(fig)



# CASOS EXIGIDOS NA ATIVIDADE

CASOS = {
    # Cada caso reúne a função, sua derivada, os chutes iniciais e a imagem de saída.
    "teste": dict(
        titulo="EXEMPLO DE TESTE: f(x) = cos(x) - x, [0, 1]",
        f=lambda x: math.cos(x) - x,
        df=lambda x: -math.sin(x) - 1,
        a=0.0, b=1.0, x0_newton=0.5, x0_sec=0.0, x1_sec=1.0,
        dominio=None, raiz_ref=0.7390851332, arquivo="grafico_teste.png"),
    "p1": dict(
        titulo="PROBLEMA 1 (Resfriamento): f(t) = 20 + 90e^(-0,15t) - 50, [0, 10]",
        f=lambda t: 20 + 90 * math.exp(-0.15 * t) - 50,
        df=lambda t: -13.5 * math.exp(-0.15 * t),
        a=0.0, b=10.0, x0_newton=0.0, x0_sec=0.0, x1_sec=10.0,
        dominio=(0.0, math.inf), raiz_ref=math.log(3) / 0.15, arquivo="grafico_problema1.png"),
    "p2": dict(
        titulo="PROBLEMA 2 (Capacitor RC): f(t) = 5(1 - e^(-t)) - 3,8, [0, 5]",
        f=lambda t: 5 * (1 - math.exp(-t)) - 3.8,
        df=lambda t: 5 * math.exp(-t),
        a=0.0, b=5.0, x0_newton=0.0, x0_sec=0.0, x1_sec=5.0,
        dominio=(0.0, math.inf), raiz_ref=-math.log(0.24), arquivo="grafico_problema2.png"),
    "p2b": dict(
        titulo="PROBLEMA 2 (Capacitor RC) - intervalo alternativo [0, 3]",
        f=lambda t: 5 * (1 - math.exp(-t)) - 3.8,
        df=lambda t: 5 * math.exp(-t),
        a=0.0, b=3.0, x0_newton=0.0, x0_sec=0.0, x1_sec=3.0,
        dominio=(0.0, math.inf), raiz_ref=-math.log(0.24), arquivo="grafico_problema2_alt.png"),
}


def executar_caso(chave, grafico=True, mostrar=True):
    """Executa um caso cadastrado, imprime a comparação e opcionalmente plota."""
    c = CASOS[chave]
    res = comparar(c["f"], c["df"], c["a"], c["b"], c["x0_newton"], c["x0_sec"],
                   c["x1_sec"], dominio=c["dominio"])
    imprimir_tabela(c["titulo"], res, c["raiz_ref"])
    print(f"Raiz de referência (analítica): {c['raiz_ref']:.10f}")
    if grafico:
        plotar(c["titulo"], c["f"], c["a"], c["b"], res, c["arquivo"], mostrar)
    return res


def demonstrar_falhas():
    """Demonstra os tratamentos de erro exigidos (sem mudança de sinal, derivada nula, etc.)."""
    print("\n" + "=" * 92)
    print("DEMONSTRAÇÃO DE FALHAS E TRATAMENTO DE EXCEÇÕES")
    print("=" * 92)
    # Usa exemplos controlados para mostrar como os métodos informam situações inválidas.
    f = lambda x: math.cos(x) - x
    df = lambda x: -math.sin(x) - 1

    print("\n[1] Ausência de mudança de sinal: cos(x)-x em [2, 3]")
    print("   ", bisseccao(f, 2, 3).mensagem)
    print("   ", falsa_posicao(f, 2, 3).mensagem)

    print("\n[2] Derivada nula: f(x) = x² - 2, f'(x) = 2x, x0 = 0")
    print("   ", newton_raphson(lambda x: x * x - 2, lambda x: 2 * x, 0.0).mensagem)

    print("\n[3] Divisão por zero na Secante: f(x) = x² - 2, x0 = -1, x1 = 1 (f(x0) = f(x1))")
    print("   ", secante(lambda x: x * x - 2, -1.0, 1.0).mensagem)

    print("\n[4] Extrapolação do domínio (Problema 2, Newton com t0 = 5, exige t >= 0)")
    c = CASOS["p2"]
    print("   ", newton_raphson(c["f"], c["df"], 5.0, dominio=(0, math.inf)).mensagem)
    print("    (Sem a verificação de domínio, o iterado cairia em t ≈ -29,6 e o exp() cresceria sem controle.)")

    print("\n[5] Divergência de Newton: f(x) = atan(x), x0 = 2 (|x_k| cresce a cada passo)")
    r = newton_raphson(math.atan, lambda x: 1 / (1 + x * x), 2.0, max_iter=20)
    print("   ", r.mensagem, "-> o método diverge por estar longe da raiz x = 0.")



# ENTRADA FLEXÍVEL DE DADOS (função digitada pelo usuário)

_NS_SEGURO = {k: getattr(math, k) for k in dir(math) if not k.startswith("_")}
# Adiciona nomes usados com frequência em português e constantes matemáticas.
_NS_SEGURO.update({"sen": math.sin, "ln": math.log, "abs": abs, "e": math.e, "pi": math.pi})


def compilar(expr, var="x"):
    """Converte texto (ex.: 'cos(x) - x', 'x**3 - 9*x + 3', 'exp(-x)*sen(x)') em função."""
    expr = expr.replace("^", "**")  # aceita a notação usual de potência
    codigo = compile(expr, "<f>", "eval")
    # A expressão só recebe a variável e os nomes permitidos em _NS_SEGURO.
    return lambda v: eval(codigo, {"__builtins__": {}}, {**_NS_SEGURO, var: v})


def ler_float(msg, padrao=None):
    """Lê um número; Enter usa o valor padrão quando ele foi fornecido."""
    while True:
        txt = input(msg).strip().replace(",", ".")
        if not txt and padrao is not None:
            return padrao
        try:
            return float(eval(txt, {"__builtins__": {}}, _NS_SEGURO))
        except Exception:
            print("  Valor inválido. Tente novamente.")


def caso_personalizado():
    """Lê a função e os parâmetros digitados e compara os métodos para esse caso."""
    global EPS_PADRAO, MAX_ITER_PADRAO
    print("\n--- Nova equação (use x como variável; funções: sin/sen, cos, exp, log/ln, sqrt, ^) ---")
    try:
        f = compilar(input("f(x) = "))
        f(1.0)
    except Exception as e:
        print(f"Expressão inválida: {e}")
        return
    deriv = input("f'(x) [ENTER = derivada numérica]: ").strip()
    if deriv:
        try:
            df = compilar(deriv); df(1.0)
        except Exception as e:
            print(f"Derivada inválida: {e}"); return
    else:
        h = 1e-7
        # Estima f'(x) pela diferença central quando o usuário não fornece a derivada.
        df = lambda x: (f(x + h) - f(x - h)) / (2 * h)   # diferença central
    a = ler_float("Início do intervalo a: ")
    b = ler_float("Fim do intervalo b: ")
    x0n = ler_float(f"x0 para Newton [ENTER = {a}]: ", a)
    x0s = ler_float(f"x0 para Secante [ENTER = {a}]: ", a)
    x1s = ler_float(f"x1 para Secante [ENTER = {b}]: ", b)
    EPS_PADRAO = ler_float(f"Precisão ε [ENTER = {EPS_PADRAO:g}]: ", EPS_PADRAO)
    MAX_ITER_PADRAO = int(ler_float(f"Máx. iterações [ENTER = {MAX_ITER_PADRAO}]: ", MAX_ITER_PADRAO))
    try:
        res = comparar(f, df, a, b, x0n, x0s, x1s, EPS_PADRAO, MAX_ITER_PADRAO)
    except Exception as e:
        print(f"Erro ao avaliar a função: {e}"); return
    imprimir_tabela("EQUAÇÃO PERSONALIZADA", res)
    if input("Plotar gráfico? (s/n): ").lower().startswith("s"):
        try:
            plotar("Equação personalizada", f, a, b, res, "grafico_personalizado.png")
        except Exception as e:
            print(f"Não foi possível plotar: {e}")



# MENU PRINCIPAL

def menu():
    """Repete o menu interativo até a pessoa escolher sair."""
    while True:
        print("\n" + "#" * 60)
        print("  SIMULADOR DE ZEROS DE FUNÇÕES - CÁLCULO NUMÉRICO")
        print("#" * 60)
        print(" 1) Exemplo de teste:  cos(x) - x")
        print(" 2) Problema 1:        Resfriamento de componente")
        print(" 3) Problema 2:        Carregamento de capacitor RC")
        print(" 4) Inserir nova equação\n 7) Problema 2 com intervalo alternativo [0, 3]")
        print(" 5) Demonstrar tratamento de falhas")
        print(" 6) Rodar tudo (1, 2, 3 e 5)")
        print(" 0) Sair")
        op = input("Opção: ").strip()
        if op == "1": executar_caso("teste")
        elif op == "2": executar_caso("p1")
        elif op == "3": executar_caso("p2")
        elif op == "4": caso_personalizado()
        elif op == "7": executar_caso("p2b")
        elif op == "5": demonstrar_falhas()
        elif op == "6":
            for k in ("teste", "p1", "p2", "p2b"): executar_caso(k)
            demonstrar_falhas()
        elif op == "0":
            print("Encerrando. Até logo!"); break
        else:
            print("Opção inválida.")
            continue
        if input("\nDeseja realizar novo cálculo? (s/n): ").strip().lower().startswith("n"):
            print("Encerrando. Até logo!"); break


if __name__ == "__main__":
    # Este bloco só roda quando o arquivo é executado diretamente, não quando é importado.
    if "--auto" in sys.argv:
        # No modo automático, cria os PNGs sem abrir janelas; --sem-grafico desativa os plots.
        graf = "--sem-grafico" not in sys.argv
        for k in ("teste", "p1", "p2", "p2b"):
            executar_caso(k, grafico=graf, mostrar=False)
        demonstrar_falhas()
    else:
        menu()
