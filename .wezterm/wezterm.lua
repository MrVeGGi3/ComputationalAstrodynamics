-- ============================================================
--  WezTerm config — ComputationalAstrodynamics
--  Layout: Editor (top-left) | Julia REPL (top-right)
--          Docker logs (bottom-left) | Claude Code (bottom-right)
-- ============================================================

local wezterm = require("wezterm")
local act     = wezterm.action
local config  = wezterm.config_builder()

-- ── Aparência ────────────────────────────────────────────────
config.color_scheme        = "Tokyo Night"
config.font                = wezterm.font("JetBrainsMono Nerd Font", { weight = "Regular" })
config.font_size           = 13.0
config.line_height         = 1.2
config.window_background_opacity = 0.96
config.text_background_opacity   = 1.0
config.enable_tab_bar      = true
config.hide_tab_bar_if_only_one_tab = false
config.use_fancy_tab_bar   = false
config.tab_bar_at_bottom   = true
config.window_padding      = { left = 8, right = 8, top = 6, bottom = 0 }

-- ── Shell padrão ─────────────────────────────────────────────
config.default_prog = { "/bin/bash", "-l" }

-- ── Variáveis de ambiente do projeto ─────────────────────────
config.set_environment_variables = {
  ASTRO_ROOT    = os.getenv("HOME") .. "/Documentos/ComputationalAstrodynamics",
  JULIA_PROJECT = os.getenv("HOME") .. "/Documentos/ComputationalAstrodynamics/julia",
  TERM          = "xterm-256color",
}

-- ── Atalhos de teclado ───────────────────────────────────────
config.keys = {
  -- Divisão de painéis
  { key = "d", mods = "CTRL|SHIFT", action = act.SplitHorizontal  { domain = "CurrentPaneDomain" } },
  { key = "e", mods = "CTRL|SHIFT", action = act.SplitVertical    { domain = "CurrentPaneDomain" } },
  -- Navegar painéis
  { key = "h", mods = "ALT",        action = act.ActivatePaneDirection "Left"  },
  { key = "l", mods = "ALT",        action = act.ActivatePaneDirection "Right" },
  { key = "k", mods = "ALT",        action = act.ActivatePaneDirection "Up"    },
  { key = "j", mods = "ALT",        action = act.ActivatePaneDirection "Down"  },
  -- Redimensionar painéis
  { key = "H", mods = "ALT|SHIFT",  action = act.AdjustPaneSize { "Left",  5 } },
  { key = "L", mods = "ALT|SHIFT",  action = act.AdjustPaneSize { "Right", 5 } },
  { key = "K", mods = "ALT|SHIFT",  action = act.AdjustPaneSize { "Up",    3 } },
  { key = "J", mods = "ALT|SHIFT",  action = act.AdjustPaneSize { "Down",  3 } },
  -- Nova aba / fechar
  { key = "t", mods = "CTRL|SHIFT", action = act.SpawnTab "CurrentPaneDomain" },
  { key = "w", mods = "CTRL|SHIFT", action = act.CloseCurrentPane { confirm = true } },
  -- Zoom no painel atual
  { key = "z", mods = "CTRL|SHIFT", action = act.TogglePaneZoomState },
  -- Abrir projeto direto
  { key = "p", mods = "CTRL|SHIFT", action = act.SpawnCommandInNewTab {
      args = { "/bin/bash", "-c", "cd $ASTRO_ROOT && exec bash" },
    }
  },
}

-- ── Layout automático ao iniciar ─────────────────────────────
-- Cria 4 painéis: editor | julia | docker | claude
wezterm.on("gui-startup", function(cmd)
  local tab, pane, window = wezterm.mux.spawn_window(cmd or {})
  local root = os.getenv("HOME") .. "/Documentos/ComputationalAstrodynamics"

  -- Painel 1 – editor / shell principal (top-left)
  pane:send_text("cd " .. root .. " && clear\n")
  tab:set_title("🛰 AstroDyn")

  -- Painel 2 – Julia REPL (top-right)
  local julia_pane = pane:split_pane {
    direction    = "Right",
    size         = 0.5,
    top_level    = false,
  }
  julia_pane:send_text(
    "cd " .. root .. "/julia && "
    .. "docker compose -f " .. root .. "/docker/docker-compose.yml "
    .. "exec julia julia --project=. -e 'using Revise; println(\"Julia pronto ✓\")'\n"
  )

  -- Painel 3 – Docker logs (bottom-left)
  local log_pane = pane:split_pane {
    direction  = "Bottom",
    size       = 0.3,
    top_level  = false,
  }
  log_pane:send_text(
    "cd " .. root .. " && "
    .. "docker compose -f docker/docker-compose.yml logs -f\n"
  )

  -- Painel 4 – Claude Code (bottom-right)
  local claude_pane = julia_pane:split_pane {
    direction  = "Bottom",
    size       = 0.3,
    top_level  = false,
  }
  claude_pane:send_text("cd " .. root .. " && claude\n")
end)

-- ── Status bar ───────────────────────────────────────────────
wezterm.on("update-right-status", function(window, _)
  local cells = {}
  local cwd   = window:active_pane():get_current_working_dir()
  if cwd then
    local path = cwd.file_path or tostring(cwd)
    path = path:gsub(os.getenv("HOME"), "~")
    table.insert(cells, "📁 " .. path)
  end
  table.insert(cells, wezterm.strftime("🕐 %H:%M"))
  window:set_right_status(wezterm.format({
    { Foreground = { AnsiColor = "Silver" } },
    { Text = table.concat(cells, "  ") },
  }))
end)

return config
