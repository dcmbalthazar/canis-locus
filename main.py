from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from genetics.engine import calculate

app = FastAPI(title="Canis Locus API")

# Servir arquivos estáticos (HTML/CSS/JS)
app.mount("/static", StaticFiles(directory="static"), name="static")


class Parent(BaseModel):
    E: tuple[str, str]
    A: tuple[str, str]
    B: tuple[str, str]
    D: tuple[str, str]
    M: tuple[str, str]
    S: tuple[str, str]


class CrossRequest(BaseModel):
    parent1: Parent
    parent2: Parent


@app.get("/")
def home():
    # Abre a tela do app
    return FileResponse("static/index.html")


@app.post("/cross")
def cross(data: CrossRequest):
    try:
        return calculate(data.parent1.model_dump(), data.parent2.model_dump())
    except ValueError as e:
        # Mensagem amigável para aparecer como erro 400 no frontend
        raise HTTPException(status_code=400, detail=str(e))
