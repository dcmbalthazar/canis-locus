from collections import defaultdict
from genetics.phenotype import interpret_phenotype


# Ordem (para impressão consistente)
LOCUS_ORDER = ["E", "A", "B", "D", "M", "S"]

# Ordem de dominância no locus A (para gerar genótipos ordenados)
A_DOM = {"Ay": 0, "at": 1, "a": 2}


def _sorted_pair(locus: str, a1: str, a2: str) -> tuple[str, str]:
    """Normaliza a ordem do par de alelos para imprimir de forma consistente."""
    if locus == "A":
        return (a1, a2) if A_DOM.get(a1, 99) <= A_DOM.get(a2, 99) else (a2, a1)
    return (a1, a2) if a1 <= a2 else (a2, a1)


def _genotype_signature(geno: dict) -> str:
    """Assinatura única para o genótipo completo (para deduplicar)."""
    parts = []
    for locus in LOCUS_ORDER:
        p = _sorted_pair(locus, geno[locus][0], geno[locus][1])
        parts.append(f"{locus}:{p[0]}/{p[1]}")
    return " | ".join(parts)


def _format_full_genotype(geno: dict) -> str:
    """Impressão completa."""
    parts = []
    for locus in LOCUS_ORDER:
        a1, a2 = _sorted_pair(locus, geno[locus][0], geno[locus][1])
        parts.append(f"{locus}:{a1}/{a2}")
    return " – ".join(parts)


def _portador_tags(geno: dict) -> list[str]:
    """
    Tags de portador para recessivos relevantes.
    """
    tags = []

    A = geno["A"]
    if "Ay" in A and ("at" in A or "a" in A):
        tags.append("Agouti(A)")

    if geno["B"].count("b") == 1:
        tags.append("marrom")

    if geno["D"].count("d") == 1:
        tags.append("diluição")

    return tags


def _validate_cross(parent1: dict, parent2: dict) -> None:
    """
    Trava anti-duplo-merle (regra do projeto):
    - Proíbe qualquer cruzamento envolvendo um duplo merle (M/M).
    - Proíbe qualquer cruzamento onde ambos os pais tenham alelo M (merle x merle),
      pois pode gerar M/M.
    Resultado: fica impossível planejar/simular filhotes duplo merle.
    """
    p1_M = parent1["M"].count("M")
    p2_M = parent2["M"].count("M")

    # 1) Se algum for M/M, bloqueia sempre
    if p1_M == 2 or p2_M == 2:
        raise ValueError(
            "Cruzamento bloqueado: um dos pais é duplo merle (M/M). "
            "O app não permite cruzar duplo merle."
        )

    # 2) Se ambos têm ao menos um M, é merle x merle (risco de M/M) -> bloqueia
    if p1_M >= 1 and p2_M >= 1:
        raise ValueError(
            "Cruzamento bloqueado: merle × merle (ambos M/_). "
            "O app bloqueia porque pode produzir filhotes duplo merle (M/M)."
        )


def _child_genotypes(parent1: dict, parent2: dict):
    """
    Gera todos os genótipos possíveis de filhotes (com probabilidades exatas).
    """
    children = [({"E": None, "A": None, "B": None, "D": None, "M": None, "S": None}, 1.0)]

    for locus in LOCUS_ORDER:
        new_children = []
        p1 = parent1[locus]
        p2 = parent2[locus]

        for child, prob in children:
            for a1 in p1:
                for a2 in p2:
                    new_child = dict(child)
                    new_child[locus] = [a1, a2]
                    new_children.append((new_child, prob * 0.25))
        children = new_children

    return children


def calculate(parent1: dict, parent2: dict):
    """
    Retorna lista de fenótipos prováveis, do mais provável -> menos provável.
    """
    # ✅ trava antes de qualquer cálculo
    _validate_cross(parent1, parent2)

    children = _child_genotypes(parent1, parent2)

    bucket_prob = defaultdict(float)
    bucket_genotypes = defaultdict(list)
    bucket_alerts = defaultdict(set)
    seen_per_bucket = defaultdict(set)

    for geno, p in children:
        phenotype, ph_alerts = interpret_phenotype(geno)

        bucket_prob[phenotype] += p

        sig = _genotype_signature(geno)
        if sig not in seen_per_bucket[phenotype]:
            seen_per_bucket[phenotype].add(sig)

            full = _format_full_genotype(geno)

            tags = _portador_tags(geno)
            if tags:
                full = f"{full} – portador:{','.join(tags)}"

            bucket_genotypes[phenotype].append(full)

        for a in ph_alerts:
            bucket_alerts[phenotype].add(a)

    results = []
    for phen, p in bucket_prob.items():
        results.append({
            "phenotype": phen,
            "probability": round(p * 100, 2),
            "genotypes": bucket_genotypes[phen],
            "alerts": sorted(list(bucket_alerts[phen])) if bucket_alerts[phen] else []
        })

    results.sort(key=lambda x: x["probability"], reverse=True)
    return results
