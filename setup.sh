#!/usr/bin/env bash
# ============================================================
#  setup.sh — ComputationalAstrodynamics
#  Idempotente: cada etapa verifica se já foi concluída.
#  Pode ser re-executado quantas vezes quiser com segurança.
#  Uso: bash setup.sh
# ============================================================

set -euo pipefail

ASTRO_ROOT="$HOME/Documentos/ComputationalAstrodynamics"
COMPOSE_FILE="$ASTRO_ROOT/docker/docker-compose.yml"
WEZTERM_CONFIG_DIR="$HOME/.config/wezterm"
WEZTERM_CONFIG_FILE="$WEZTERM_CONFIG_DIR/wezterm.lua"
WEZTERM_SOURCE="$ASTRO_ROOT/.wezterm/wezterm.lua"

BOLD='\033[1m'; CYAN='\033[0;36m'; GREEN='\033[0;32m'
YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'

info()    { echo -e "${CYAN}[INFO]${NC}  $*"; }
success() { echo -e "${GREEN}[OK]${NC}    $*"; }
skip()    { echo -e "${GREEN}[SKIP]${NC}  $* (já concluído)"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}  $*"; }
fail()    { echo -e "${RED}[FALHA]${NC} $*"; FAILURES+=("$*"); }

FAILURES=()

echo -e "\n${BOLD}🛰  ComputationalAstrodynamics — Setup${NC}"
echo -e "   Re-executável com segurança — etapas já concluídas são puladas.\n"

# ════════════════════════════════════════════════════════════
# ETAPA 1 — Estrutura de pastas
# ════════════════════════════════════════════════════════════
echo -e "${BOLD}[1/8] Estrutura de pastas${NC}"
DIRS=(
    "$ASTRO_ROOT/.wezterm"
    "$ASTRO_ROOT/docker"
    "$ASTRO_ROOT/julia/src"
    "$ASTRO_ROOT/julia/scripts"
    "$ASTRO_ROOT/julia/notebooks"
    "$ASTRO_ROOT/julia/data/input"
    "$ASTRO_ROOT/julia/data/output"
    "$ASTRO_ROOT/claude"
    "$ASTRO_ROOT/.config"
)
ALL_EXIST=true
for d in "${DIRS[@]}"; do [ -d "$d" ] || ALL_EXIST=false; done

if $ALL_EXIST; then
    skip "Todas as pastas existem"
else
    mkdir -p "${DIRS[@]}"
    success "Pastas criadas"
fi

# ════════════════════════════════════════════════════════════
# ETAPA 2 — docker-compose.yml presente e válido
# ════════════════════════════════════════════════════════════
echo -e "\n${BOLD}[2/8] docker-compose.yml${NC}"
if [ ! -f "$COMPOSE_FILE" ]; then
    fail "docker-compose.yml não encontrado em $COMPOSE_FILE"
    warn "Verifique se os arquivos do projeto foram copiados para $ASTRO_ROOT"
else
    if command -v python3 &>/dev/null; then
        if python3 -c "import yaml; yaml.safe_load(open('$COMPOSE_FILE'))" 2>/dev/null; then
            success "docker-compose.yml encontrado e sintaxe YAML válida"
        else
            fail "docker-compose.yml tem erro de sintaxe YAML"
        fi
    else
        success "docker-compose.yml encontrado"
    fi
    if grep -q "julia:" "$COMPOSE_FILE"; then
        success "Serviço 'julia' declarado no compose"
    else
        fail "Serviço 'julia' não encontrado no docker-compose.yml"
    fi
fi

# ════════════════════════════════════════════════════════════
# ETAPA 3 — WezTerm
# ════════════════════════════════════════════════════════════
echo -e "\n${BOLD}[3/8] WezTerm${NC}"
if [ ! -f "$WEZTERM_SOURCE" ]; then
    fail "Arquivo fonte do WezTerm não encontrado: $WEZTERM_SOURCE"
else
    mkdir -p "$WEZTERM_CONFIG_DIR"
    if [ -L "$WEZTERM_CONFIG_FILE" ] && [ "$(readlink -f "$WEZTERM_CONFIG_FILE")" = "$WEZTERM_SOURCE" ]; then
        skip "Symlink já aponta para $WEZTERM_SOURCE"
    else
        if [ -f "$WEZTERM_CONFIG_FILE" ] && [ ! -L "$WEZTERM_CONFIG_FILE" ]; then
            cp "$WEZTERM_CONFIG_FILE" "${WEZTERM_CONFIG_FILE}.bak"
            warn "Backup do wezterm.lua anterior → ${WEZTERM_CONFIG_FILE}.bak"
        fi
        ln -sf "$WEZTERM_SOURCE" "$WEZTERM_CONFIG_FILE"
        success "Symlink criado: $WEZTERM_CONFIG_FILE → $WEZTERM_SOURCE"
    fi
fi

# ════════════════════════════════════════════════════════════
# ETAPA 4 — Docker daemon
# ════════════════════════════════════════════════════════════
echo -e "\n${BOLD}[4/8] Docker${NC}"
DOCKER_OK=false
if ! command -v docker &>/dev/null; then
    fail "Docker não instalado"
    warn "Instale em: https://docs.docker.com/engine/install/ubuntu/"
    warn "Após instalar: sudo usermod -aG docker \$USER && newgrp docker"
else
    DOCKER_VER=$(docker --version | awk '{print $3}' | tr -d ',')
    if docker info &>/dev/null 2>&1; then
        success "Docker $DOCKER_VER ativo e acessível sem sudo"
        DOCKER_OK=true
    else
        fail "Docker instalado ($DOCKER_VER) mas daemon inacessível"
        warn "• Daemon parado?  → sudo systemctl start docker"
        warn "• Sem permissão?  → sudo usermod -aG docker \$USER  (requer logout/login)"
    fi
fi

# ════════════════════════════════════════════════════════════
# ETAPA 5 — Docker Compose plugin + validação do compose file
# ════════════════════════════════════════════════════════════
echo -e "\n${BOLD}[5/8] Docker Compose${NC}"
COMPOSE_OK=false
if ! $DOCKER_OK; then
    warn "Pulando (Docker com problema na etapa anterior)"
elif ! docker compose version &>/dev/null 2>&1; then
    fail "Docker Compose plugin não encontrado"
    warn "Instale: sudo apt install docker-compose-plugin"
else
    COMPOSE_VER=$(docker compose version --short 2>/dev/null || \
                  docker compose version | grep -oP '\d+\.\d+\.\d+' | head -1)
    if docker compose -f "$COMPOSE_FILE" config --quiet 2>/dev/null; then
        success "Docker Compose $COMPOSE_VER — compose file validado com sucesso"
        COMPOSE_OK=true
    else
        fail "docker compose config falhou — erro no $COMPOSE_FILE"
        info "Diagnóstico:"
        docker compose -f "$COMPOSE_FILE" config 2>&1 | head -20 || true
    fi
fi

# ════════════════════════════════════════════════════════════
# ETAPA 6 — Imagem Docker
# ════════════════════════════════════════════════════════════
echo -e "\n${BOLD}[6/8] Imagem Docker${NC}"
if ! $COMPOSE_OK; then
    warn "Pulando build (Docker Compose com problema na etapa anterior)"
elif docker image inspect astrodynamics-julia:latest &>/dev/null 2>&1; then
    IMG_DATE=$(docker image inspect astrodynamics-julia:latest \
               --format '{{.Created}}' | cut -c1-10)
    skip "Imagem astrodynamics-julia:latest existe (criada $IMG_DATE)"
    info "Para reconstruir: docker compose -f $COMPOSE_FILE build --no-cache"
else
    echo ""
    read -r -p "  Imagem não encontrada. Construir agora? (~5 min) [s/N] " BUILD_NOW
    if [[ "$BUILD_NOW" =~ ^[Ss]$ ]]; then
        info "Construindo imagem..."
        if cd "$ASTRO_ROOT" && docker compose -f "$COMPOSE_FILE" build 2>&1; then
            success "Imagem construída com sucesso"
        else
            fail "Falha no docker build — veja erros acima"
        fi
    else
        warn "Build pulado. Execute depois:"
        echo "    docker compose -f $COMPOSE_FILE build"
    fi
fi

# ════════════════════════════════════════════════════════════
# ETAPA 7 — Claude Code
# ════════════════════════════════════════════════════════════
echo -e "\n${BOLD}[7/8] Claude Code${NC}"
export PATH="$HOME/.local/bin:$PATH"
if command -v claude &>/dev/null; then
    CLAUDE_VER=$(claude --version 2>/dev/null || echo "versão desconhecida")
    skip "Claude Code já instalado: $CLAUDE_VER"
    info "Para atualizar: claude update"
else
    info "Baixando e instalando Claude Code (instalador nativo)..."
    if curl -fsSL https://claude.ai/install.sh | bash; then
        export PATH="$HOME/.local/bin:$PATH"
        if command -v claude &>/dev/null; then
            success "Claude Code instalado com sucesso"
        else
            fail "'claude' não encontrado no PATH após instalação"
            warn "Adicione ao ~/.bashrc:  export PATH=\"\$HOME/.local/bin:\$PATH\""
        fi
    else
        fail "Falha no download/instalação do Claude Code"
        warn "Tente manualmente: curl -fsSL https://claude.ai/install.sh | bash"
    fi
fi

# ════════════════════════════════════════════════════════════
# ETAPA 8 — Variáveis de ambiente
# ════════════════════════════════════════════════════════════
echo -e "\n${BOLD}[8/8] Variáveis de ambiente${NC}"
SHELL_RC=""
[ -f "$HOME/.zshrc" ]  && SHELL_RC="$HOME/.zshrc"
[ -f "$HOME/.bashrc" ] && SHELL_RC="$HOME/.bashrc"

MARKER="# ComputationalAstrodynamics"
if [ -z "$SHELL_RC" ]; then
    warn "Nenhum ~/.bashrc ou ~/.zshrc encontrado"
elif grep -q "$MARKER" "$SHELL_RC"; then
    skip "Variáveis já presentes em $SHELL_RC"
else
    cat >> "$SHELL_RC" <<EOF

$MARKER
export ASTRO_ROOT="$ASTRO_ROOT"
export JULIA_PROJECT="\$ASTRO_ROOT/julia"
export JULIA_NUM_THREADS="auto"
export JULIA_REVISE_INCLUDE="1"
export PATH="\$HOME/.local/bin:\$PATH"
EOF
    success "Variáveis adicionadas a $SHELL_RC"
fi

# ════════════════════════════════════════════════════════════
# RESUMO FINAL
# ════════════════════════════════════════════════════════════
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

if [ ${#FAILURES[@]} -eq 0 ]; then
    echo -e "${GREEN}${BOLD}  ✓ Todas as etapas concluídas com sucesso!${NC}"
else
    echo -e "${RED}${BOLD}  ✗ ${#FAILURES[@]} problema(s) encontrado(s):${NC}"
    for f in "${FAILURES[@]}"; do
        echo -e "    ${RED}•${NC} $f"
    done
    echo ""
    echo -e "  Corrija os problemas e rode ${BOLD}bash setup.sh${NC} novamente."
    echo -e "  As etapas já concluídas serão puladas automaticamente."
fi

echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
echo "  Projeto:  $ASTRO_ROOT"
echo "  Compose:  $COMPOSE_FILE"
echo "  WezTerm:  $WEZTERM_CONFIG_FILE"
echo ""

if [ ${#FAILURES[@]} -eq 0 ]; then
    echo -e "${BOLD}  Próximos passos:${NC}"
    echo "  1. source $SHELL_RC          (ou abra novo terminal)"
    echo "  2. cd \$ASTRO_ROOT && claude  (login Pro/Max no browser)"
    echo "  3. Abra o WezTerm            (4 painéis criados automaticamente)"
fi
echo ""
